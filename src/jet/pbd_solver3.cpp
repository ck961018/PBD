// Copyright (c) 2018 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#include <pch.h>

#include <jet/array_utils.h>
#include <jet/constant_vector_field3.h>
#include <jet/parallel.h>
#include <jet/pbd_solver3.h>
#include <jet/timer.h>

#include <algorithm>

namespace jet {

PbdSolver3::PbdSolver3() : PbdSolver3(1e-3) {}

PbdSolver3::PbdSolver3(double mass) {
  _pbdData = std::make_shared<PbdData3>();
  _pbdData->setMass(mass);
  _wind = std::make_shared<ConstantVectorField3>(Vector3D());
}

PbdSolver3::~PbdSolver3() {}

double PbdSolver3::dragCoefficient() const { return _dragCoefficient; }

void PbdSolver3::setDragCoefficient(double newDragCoefficient) {
  _dragCoefficient = std::max(newDragCoefficient, 0.0);
}

double PbdSolver3::restitutionCoefficient() const { return _restitutionCoefficient; }

void PbdSolver3::setRestitutionCoefficient(double newRestitutionCoefficient) {
  _restitutionCoefficient = clamp(newRestitutionCoefficient, 0.0, 1.0);
}

const Vector3D& PbdSolver3::gravity() const { return _gravity; }

void PbdSolver3::setGravity(const Vector3D& newGravity) { _gravity = newGravity; }

const double& PbdSolver3::edgeCompliance() const { return _edgeCompliance; }

void PbdSolver3::setEdgeCompliance(const double& newEdgeCompliance) {
  _edgeCompliance = newEdgeCompliance;
}

const double& PbdSolver3::volumeCompliance() const { return _volumeCompliance; }

void PbdSolver3::setVolumeCompliance(const double& newVolumeCompliance) {
  _volumeCompliance = newVolumeCompliance;
}

const PbdData3Ptr& PbdSolver3::pbdData() const { return _pbdData; }

const Collider3Ptr& PbdSolver3::collider() const { return _collider; }

void PbdSolver3::setCollider(const Collider3Ptr& newCollider) { _collider = newCollider; }

const VectorField3Ptr& PbdSolver3::wind() const { return _wind; }

void PbdSolver3::setWind(const VectorField3Ptr& newWind) { _wind = newWind; }

void PbdSolver3::onInitialize() {
  // When initializing the solver, update the collider and emitter state as
  // well since they also affects the initial condition of the simulation.
  Timer timer;
  updateCollider(0.0);
  JET_INFO << "Update collider took " << timer.durationInSeconds() << " seconds";
}

void PbdSolver3::onAdvanceTimeStep(double timeStepInSeconds) {
  beginAdvanceTimeStep(timeStepInSeconds);

  Timer timer;
  accumulateForces(timeStepInSeconds);
  JET_INFO << "Accumulating forces took " << timer.durationInSeconds() << " seconds";

  timer.reset();
  timeIntegration(timeStepInSeconds);
  JET_INFO << "Time integration took " << timer.durationInSeconds() << " seconds";

  timer.reset();
  resolveCollision();
  JET_INFO << "Resolving collision took " << timer.durationInSeconds() << " seconds";

  timer.reset();
  resolveConstraints(timeStepInSeconds);
  JET_INFO << "Resolving constraints took " << timer.durationInSeconds() << " seconds";

  endAdvanceTimeStep(timeStepInSeconds);
}

void PbdSolver3::accumulateForces(double timeStepInSeconds) {
  UNUSED_VARIABLE(timeStepInSeconds);

  // Add external forces
  accumulateExternalForces();
}

void PbdSolver3::beginAdvanceTimeStep(double timeStepInSeconds) {
  // Clear forces
  auto forces = _pbdData->forces();
  setRange1(forces.size(), Vector3D(), &forces);

  // Update collider
  Timer timer;
  updateCollider(timeStepInSeconds);
  JET_INFO << "Update collider took " << timer.durationInSeconds() << " seconds";

  // Allocate buffers
  size_t n = _pbdData->numberOfParticles();
  _newPositions.resize(n);
  _newVelocities.resize(n);

  onBeginAdvanceTimeStep(timeStepInSeconds);
}

void PbdSolver3::endAdvanceTimeStep(double timeStepInSeconds) {
  // Update data
  size_t n = _pbdData->numberOfParticles();
  auto positions = _pbdData->positions();
  auto velocities = _pbdData->velocities();

  parallelFor(kZeroSize, n, [&](size_t i) {
    // post solve
    velocities[i] = (_newPositions[i] - positions[i]) / timeStepInSeconds;
    positions[i] = _newPositions[i];
  });

  onEndAdvanceTimeStep(timeStepInSeconds);
}

void PbdSolver3::onBeginAdvanceTimeStep(double timeStepInSeconds) {
  UNUSED_VARIABLE(timeStepInSeconds);
}

void PbdSolver3::onEndAdvanceTimeStep(double timeStepInSeconds) {
  UNUSED_VARIABLE(timeStepInSeconds);
}

void PbdSolver3::resolveCollision() {
  resolveCollision(_newPositions.accessor(), _newVelocities.accessor());
}

void PbdSolver3::resolveCollision(ArrayAccessor1<Vector3D> newPositions,
                                  ArrayAccessor1<Vector3D> newVelocities) {
  if (_collider != nullptr) {
    size_t numberOfParticles = _pbdData->numberOfParticles();

    parallelFor(kZeroSize, numberOfParticles, [&](size_t i) {
      _collider->resolveCollision(0, _restitutionCoefficient, &newPositions[i], &newVelocities[i]);
    });
  }
}

void PbdSolver3::resolveConstraints(double timeStepInSeconds) {
  auto restLen = _pbdData->restLen();
  auto restVol = _pbdData->restVol();
  auto invMass = _pbdData->invMass();
  // solve edges
  double alpha = _edgeCompliance / timeStepInSeconds / timeStepInSeconds;
  auto positionIdxOfEdge = _pbdData->positionIdxOfEdges();

  positionIdxOfEdge.parallelForEachIndex([&](size_t i) {
    auto& _idx = positionIdxOfEdge[i];
    size_t idx[2];
    idx[0] = std::get<0>(_idx);
    idx[1] = std::get<1>(_idx);

    auto& p0 = _newPositions[idx[0]];
    auto& p1 = _newPositions[idx[1]];

    Vector3D edgeGrads = p0 - p1;
    double len = edgeGrads.length();
    edgeGrads /= len;
    double C = len - restLen[i];

    double w = invMass[idx[0]] + invMass[idx[1]];

    double s = -C / (w + alpha);

    _newPositions[idx[0]] += edgeGrads * (+s * invMass[idx[0]]);
    _newPositions[idx[1]] += edgeGrads * (-s * invMass[idx[1]]);
  });

  // solve volumes
  alpha = _volumeCompliance / timeStepInSeconds / timeStepInSeconds;
  auto positionIdxOfTet = _pbdData->positionIdxOfTets();
  positionIdxOfTet.parallelForEachIndex([&](size_t i) {
    auto& _idx = positionIdxOfTet[i];
    size_t idx[4];
    idx[0] = std::get<0>(_idx);
    idx[1] = std::get<1>(_idx);
    idx[2] = std::get<2>(_idx);
    idx[3] = std::get<3>(_idx);

    Vector3D p[4];
    for (size_t j = 0; j < 4; ++j) {
      p[j] = _newPositions[idx[j]];
    }

    Vector3D volGrads[4];
    volGrads[0] = (p[3] - p[1]).cross(p[2] - p[1]);
    volGrads[1] = (p[2] - p[0]).cross(p[3] - p[0]);
    volGrads[2] = (p[3] - p[0]).cross(p[1] - p[0]);
    volGrads[3] = (p[1] - p[0]).cross(p[2] - p[0]);

    double w = 0;
    for (size_t j = 0; j < 4; ++j) {
      w += invMass[idx[j]] * (volGrads[j].lengthSquared());
    }

    if (fabs(w) < 0.000001) {
      return;
    }

    double vol = computeVolume(p[0], p[1], p[2], p[3]);
    double C = (vol - restVol[i]) * 6.0;
    double s = -C / (w + alpha);

    for (size_t j = 0; j < 4; ++j) {
      _newPositions[idx[j]] += volGrads[j] * s * invMass[idx[j]];
    }
  });
}

void PbdSolver3::setPbdData(const PbdData3Ptr& newParticles) { _pbdData = newParticles; }

void PbdSolver3::accumulateExternalForces() {
  size_t n = _pbdData->numberOfParticles();
  auto forces = _pbdData->forces();
  auto velocities = _pbdData->velocities();
  auto positions = _pbdData->positions();
  const double mass = _pbdData->mass();

  parallelFor(kZeroSize, n, [&](size_t i) {
    // Gravity
    Vector3D force = mass * _gravity;

    // Wind forces
    Vector3D relativeVel = velocities[i] - _wind->sample(positions[i]);
    force += -_dragCoefficient * relativeVel;

    forces[i] += force;
  });
}

void PbdSolver3::timeIntegration(double timeStepInSeconds) {
  size_t n = _pbdData->numberOfParticles();
  auto forces = _pbdData->forces();
  auto velocities = _pbdData->velocities();
  auto positions = _pbdData->positions();
  const double mass = _pbdData->mass();

  parallelFor(kZeroSize, n, [&](size_t i) {
    // Integrate velocity first
    Vector3D& newVelocity = _newVelocities[i];
    newVelocity = velocities[i] + timeStepInSeconds * forces[i] / mass;

    // Integrate position.
    Vector3D& newPosition = _newPositions[i];
    newPosition = positions[i] + timeStepInSeconds * newVelocity;
  });
}

void PbdSolver3::updateCollider(double timeStepInSeconds) {
  if (_collider != nullptr) {
    _collider->update(currentTimeInSeconds(), timeStepInSeconds);
  }
}

double PbdSolver3::computeVolume(Vector3D A, Vector3D B, Vector3D C, Vector3D D) {
  auto tmp = (B - A).cross(C - A);
  return tmp.dot(D - A) / 6.0;
}

PbdSolver3::Builder PbdSolver3::builder() { return Builder(); }

PbdSolver3 PbdSolver3::Builder::build() const { return PbdSolver3(_mass); }

PbdSolver3Ptr PbdSolver3::Builder::makeShared() const {
  return std::shared_ptr<PbdSolver3>(new PbdSolver3(_mass), [](PbdSolver3* obj) { delete obj; });
}

}  // namespace jet
