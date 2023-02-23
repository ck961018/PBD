#ifndef INCLUDE_JET_PBD_SOLVER3_H_
#define INCLUDE_JET_PBD_SOLVER3_H_

#include <jet/collider3.h>
#include <jet/constants.h>
#include <jet/vector_field3.h>
#include <jet/pbd_data3.h>
#include <jet/physics_animation.h>

namespace jet {

class PbdSolver3 : public PhysicsAnimation {
 public:
  class Builder;

  //! Constructs an empty solver.
  PbdSolver3();

  //! Constructs a solver with particle parameters.
  PbdSolver3(double mass);

  //! Destructor.
  virtual ~PbdSolver3();

  //! Returns the drag coefficient.
  double dragCoefficient() const;

  //!
  //! \brief      Sets the drag coefficient.
  //!
  //! The drag coefficient controls the amount of air-drag. The coefficient
  //! should be a positive number and 0 means no drag force.
  //!
  //! \param[in]  newDragCoefficient The new drag coefficient.
  //!
  void setDragCoefficient(double newDragCoefficient);

  //! Sets the restitution coefficient.
  double restitutionCoefficient() const;

  //!
  //! \brief      Sets the restitution coefficient.
  //!
  //! The restitution coefficient controls the bouncy-ness of a particle when
  //! it hits a collider surface. The range of the coefficient should be 0 to
  //! 1 -- 0 means no bounce back and 1 means perfect reflection.
  //!
  //! \param[in]  newRestitutionCoefficient The new restitution coefficient.
  //!
  void setRestitutionCoefficient(double newRestitutionCoefficient);

  //! Returns the gravity.
  const Vector3D& gravity() const;

  //! Sets the gravity.
  void setGravity(const Vector3D& newGravity);

  const double& edgeCompliance() const;

  void setEdgeCompliance(const double& newEdgeCompliance);

  const double& volumeCompliance() const;

  void setVolumeCompliance(const double& newVolumeCompliance);

  //!
  //! \brief      Returns the particle system data.
  //!
  //! This function returns the particle system data. The data is created when
  //! this solver is constructed and also owned by the solver.
  //!
  //! \return     The particle system data.
  //!
  const PbdData3Ptr& pbdData() const;

  //! Returns the collider.
  const Collider3Ptr& collider() const;

  //! Sets the collider.
  void setCollider(const Collider3Ptr& newCollider);

  //! Returns the wind field.
  const VectorField3Ptr& wind() const;

  //!
  //! \brief      Sets the wind.
  //!
  //! Wind can be applied to the particle system by setting a vector field to
  //! the solver.
  //!
  //! \param[in]  newWind The new wind.
  //!
  void setWind(const VectorField3Ptr& newWind);

  //! Returns builder fox PbdSolver3.
  static Builder builder();

 protected:
  //! Initializes the simulator.
  void onInitialize() override;

  //! Called to advane a single time-step.
  void onAdvanceTimeStep(double timeStepInSeconds) override;

  //! Accumulates forces applied to the particles.
  virtual void accumulateForces(double timeStepInSeconds);

  //! Called when a time-step is about to begin.
  virtual void onBeginAdvanceTimeStep(double timeStepInSeconds);

  //! Called after a time-step is completed.
  virtual void onEndAdvanceTimeStep(double timeStepInSeconds);

  //! Resolves any collisions occured by the particles.
  void resolveCollision();

  //! Resolves any collisions occured by the particles where the particle
  //! state is given by the position and velocity arrays.
  void resolveCollision(ArrayAccessor1<Vector3D> newPositions,
                        ArrayAccessor1<Vector3D> newVelocities);

  void resolveConstraints(double timeStepInSeconds);

  //! Assign a new particle system data.
  void setPbdData(const PbdData3Ptr& newParticles);

 private:
  double _dragCoefficient = 1e-4;
  double _restitutionCoefficient = 0.0;
  double _edgeCompliance = 100.0;
  double _volumeCompliance = 0.0;
  Vector3D _gravity = Vector3D(0.0, kGravity, 0.0);

  PbdData3Ptr _pbdData;
  PbdData3::VectorData _newPositions;
  PbdData3::VectorData _newVelocities;
  Collider3Ptr _collider;
  VectorField3Ptr _wind;

  void beginAdvanceTimeStep(double timeStepInSeconds);

  void endAdvanceTimeStep(double timeStepInSeconds);

  void accumulateExternalForces();

  void timeIntegration(double timeStepInSeconds);

  void updateCollider(double timeStepInSeconds);

  double computeVolume(Vector3D A, Vector3D B, Vector3D C, Vector3D D);
};

//! Shared pointer type for the PbdSolver3.
typedef std::shared_ptr<PbdSolver3> PbdSolver3Ptr;

//!
//! \brief Base class for particle-based solver builder.
//!
template <typename DerivedBuilder>
class PbdSolverBuilderBase3 {
 public:
  //! Returns builder with mass per particle.
  DerivedBuilder& withMass(double mass);

 protected:
  double _mass = 1;
};

template <typename T>
T& PbdSolverBuilderBase3<T>::withMass(double mass) {
  _mass = mass;
  return static_cast<T&>(*this);
}

//!
//! \brief Front-end to create PbdSolver3 objects step by step.
//!
class PbdSolver3::Builder final : public PbdSolverBuilderBase3<PbdSolver3::Builder> {
 public:
  //! Builds PbdSolver3.
  PbdSolver3 build() const;

  //! Builds shared pointer of PbdSolver3 instance.
  PbdSolver3Ptr makeShared() const;
};

}  // namespace jet

#endif  // INCLUDE_JET_PARTICLE_SYSTEM_SOLVER3_H_
