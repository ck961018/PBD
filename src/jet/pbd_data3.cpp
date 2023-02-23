// Copyright (c) 2018 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#ifdef _MSC_VER
#pragma warning(disable : 4244)
#endif

#include <pch.h>

#include <fbs_helpers.h>
// #include <generated/pbd_system_data3_generated.h>

#include <jet/parallel.h>
#include <jet/pbd_data3.h>
#include <jet/timer.h>

#include <algorithm>
#include <vector>

using namespace jet;

PbdData3::PbdData3() : PbdData3(0) {}

PbdData3::PbdData3(size_t numberOfParticles) {
  _positionIdx = addVectorData();
  _velocityIdx = addVectorData();
  _forceIdx = addVectorData();

  _invMassIdx = addScalarData();
  _restLenIdx = addScalarData();
  _restVolIdx = addScalarData();

  resize(numberOfParticles);
}

PbdData3::PbdData3(const PbdData3& other) { set(other); }

PbdData3::~PbdData3() {}

void PbdData3::resize(size_t numberOfParticles) {
  _numberOfParticles = numberOfParticles;

  _scalarDataList[_invMassIdx].resize(numberOfParticles, 0.0);

  for (auto& attr : _vectorDataList) {
    attr.resize(numberOfParticles, Vector3D());
  }
}

size_t PbdData3::numberOfParticles() const { return _numberOfParticles; }

size_t PbdData3::addScalarData(double initialVal) {
  size_t attrIdx = _scalarDataList.size();
  _scalarDataList.emplace_back(numberOfParticles(), initialVal);
  return attrIdx;
}

size_t PbdData3::addVectorData(const Vector3D& initialVal) {
  size_t attrIdx = _vectorDataList.size();
  _vectorDataList.emplace_back(numberOfParticles(), initialVal);
  return attrIdx;
}

double PbdData3::mass() const { return _mass; }

void PbdData3::setMass(double newMass) { _mass = std::max(newMass, 0.0); }

ConstArrayAccessor1<Vector3D> PbdData3::positions() const { return vectorDataAt(_positionIdx); }

ArrayAccessor1<Vector3D> PbdData3::positions() { return vectorDataAt(_positionIdx); }

ConstArrayAccessor1<Vector3D> PbdData3::velocities() const { return vectorDataAt(_velocityIdx); }

ArrayAccessor1<Vector3D> PbdData3::velocities() { return vectorDataAt(_velocityIdx); }

ConstArrayAccessor1<Vector3D> PbdData3::forces() const { return vectorDataAt(_forceIdx); }

ArrayAccessor1<Vector3D> PbdData3::forces() { return vectorDataAt(_forceIdx); }

ConstArrayAccessor1<double> PbdData3::restLen() const { return scalarDataAt(_restLenIdx); }

ArrayAccessor1<double> PbdData3::restLen() { return scalarDataAt(_restLenIdx); }

ConstArrayAccessor1<double> PbdData3::restVol() const { return scalarDataAt(_restVolIdx); }

ArrayAccessor1<double> PbdData3::restVol() { return scalarDataAt(_restVolIdx); }

ConstArrayAccessor1<double> PbdData3::invMass() const { return scalarDataAt(_invMassIdx); }

ArrayAccessor1<double> PbdData3::invMass() { return scalarDataAt(_invMassIdx); }

ConstArrayAccessor1<double> PbdData3::scalarDataAt(size_t idx) const {
  return _scalarDataList[idx].constAccessor();
}

ArrayAccessor1<double> PbdData3::scalarDataAt(size_t idx) {
  return _scalarDataList[idx].accessor();
}

ConstArrayAccessor1<Vector3D> PbdData3::vectorDataAt(size_t idx) const {
  return _vectorDataList[idx].constAccessor();
}

ArrayAccessor1<Vector3D> PbdData3::vectorDataAt(size_t idx) {
  return _vectorDataList[idx].accessor();
}

IdxOfEdge jet::PbdData3::positionIdxOfEdge(size_t idx) { return _positionIdxOfEdges[idx]; }

ArrayAccessor1<IdxOfEdge> PbdData3::positionIdxOfEdges() { return _positionIdxOfEdges.accessor(); }

ConstArrayAccessor1<IdxOfEdge> PbdData3::positionIdxOfEdges() const {
  return _positionIdxOfEdges.constAccessor();
}

IdxOfTet jet::PbdData3::positionIdxOfTet(size_t idx) { return _positionIdxOfTets[idx]; }

ArrayAccessor1<IdxOfTet> PbdData3::positionIdxOfTets() { return _positionIdxOfTets.accessor(); }

ConstArrayAccessor1<IdxOfTet> PbdData3::positionIdxOfTets() const {
  return _positionIdxOfTets.constAccessor();
}

double jet::PbdData3::lengthOfEdge(size_t idx) {
  auto& pos = positions();
  auto& positionIdx = positionIdxOfEdge(idx);
  return (pos[std::get<0>(positionIdx)] - pos[std::get<1>(positionIdx)]).length();
}

double jet::PbdData3::volumeOfTet(size_t idx) {
  auto& pos = positions();
  auto& positionIdx = positionIdxOfTet(idx);
  auto& A = pos[std::get<0>(positionIdx)];
  auto& B = pos[std::get<1>(positionIdx)];
  auto& C = pos[std::get<2>(positionIdx)];
  auto& D = pos[std::get<3>(positionIdx)];
  double res = ((B - A).cross(C - A)).dot(D - A);
  return res / 6.0;
}

void PbdData3::addParticle(const Vector3D& newPosition, const Vector3D& newVelocity,
                           const Vector3D& newForce) {
  Array1<Vector3D> newPositions = {newPosition};
  Array1<Vector3D> newVelocities = {newVelocity};
  Array1<Vector3D> newForces = {newForce};

  addParticles(newPositions.constAccessor(), newVelocities.constAccessor(),
               newForces.constAccessor());
}

void PbdData3::addParticles(const ConstArrayAccessor1<Vector3D>& newPositions,
                            const ConstArrayAccessor1<Vector3D>& newVelocities,
                            const ConstArrayAccessor1<Vector3D>& newForces) {
  JET_THROW_INVALID_ARG_IF(newVelocities.size() > 0 && newVelocities.size() != newPositions.size());
  JET_THROW_INVALID_ARG_IF(newForces.size() > 0 && newForces.size() != newPositions.size());

  size_t oldNumberOfParticles = numberOfParticles();
  size_t numberOfParticles = oldNumberOfParticles + newPositions.size();

  resize(numberOfParticles);

  auto pos = positions();
  auto vel = velocities();
  auto frc = forces();

  parallelFor(kZeroSize, newPositions.size(),
              [&](size_t i) { pos[i + oldNumberOfParticles] = newPositions[i]; });

  if (newVelocities.size() > 0) {
    parallelFor(kZeroSize, newPositions.size(),
                [&](size_t i) { vel[i + oldNumberOfParticles] = newVelocities[i]; });
  }

  if (newForces.size() > 0) {
    parallelFor(kZeroSize, newPositions.size(),
                [&](size_t i) { frc[i + oldNumberOfParticles] = newForces[i]; });
  }
}

void jet::PbdData3::addEdge(const IdxOfEdge& positionIdxOfEdge) {
  auto edgeIdx = _positionIdxOfEdges.size();
  _positionIdxOfEdges.append(positionIdxOfEdge);
  _scalarDataList[_restLenIdx].append(lengthOfEdge(edgeIdx));
}

void jet::PbdData3::addTet(const IdxOfTet& positionIdxOfTet) {
  auto tetIdx = _positionIdxOfTets.size();
  _positionIdxOfTets.append(positionIdxOfTet);
  _scalarDataList[_restVolIdx].append(volumeOfTet(tetIdx));

  auto& _idx = positionIdxOfTet;
  size_t idx[4];
  idx[0] = std::get<0>(_idx);
  idx[1] = std::get<1>(_idx);
  idx[2] = std::get<2>(_idx);
  idx[3] = std::get<3>(_idx);

  for (size_t i = 0; i < 4; i++) {
    double pInvMass = _mass / (_scalarDataList[_restVolIdx][tetIdx] / 4.0);
    _scalarDataList[_invMassIdx][idx[i]] += pInvMass;
  }
}

void PbdData3::serialize(std::vector<uint8_t>* buffer) const {
  // flatbuffers::FlatBufferBuilder builder(1024);
  // flatbuffers::Offset<fbs::PbdData3> fbsPbdData3;

  // serializedPbdSystemData(&builder, &fbsPbdData3);

  // builder.Finish(fbsPbdData3);

  // uint8_t* buf = builder.GetBufferPointer();
  // size_t size = builder.GetSize();

  // buffer->resize(size);
  // memcpy(buffer->data(), buf, size);
}

void PbdData3::deserialize(const std::vector<uint8_t>& buffer) {
  // auto fbsPbdData3 = fbs::GetPbdData3(buffer.data());
  // deserializedPbdSystemData(fbsPbdData3);
}

void PbdData3::set(const PbdData3& other) {
  _mass = other._mass;
  _positionIdx = other._positionIdx;
  _velocityIdx = other._velocityIdx;
  _forceIdx = other._forceIdx;
  _numberOfParticles = other._numberOfParticles;

  for (auto& attr : other._scalarDataList) {
    _scalarDataList.emplace_back(attr);
  }

  for (auto& attr : other._vectorDataList) {
    _vectorDataList.emplace_back(attr);
  }
}

PbdData3& PbdData3::operator=(const PbdData3& other) {
  set(other);
  return *this;
}

// void PbdData3::serializedPbdSystemData(
//     flatbuffers::FlatBufferBuilder* builder,
//     flatbuffers::Offset<fbs::PbdData3>* fbsPbdData3) const {
//   // Copy data
//   std::vector<flatbuffers::Offset<fbs::ScalardPbdData3>> scalarDataList;
//   for (const auto& scalarData : _scalarDataList) {
//     auto fbsScalarData = fbs::CreateScalardPbdData3(
//         *builder, builder->CreateVector(scalarData.data(), scalarData.size()));
//     scalarDataList.push_back(fbsScalarData);
//   }
//   auto fbsScalarDataList = builder->CreateVector(scalarDataList);

//   std::vector<flatbuffers::Offset<fbs::VectordPbdData3>> vectorDataList;
//   for (const auto& vectorData : _vectorDataList) {
//     std::vector<fbs::Vector3D> newVectorData;
//     for (const auto& v : vectorData) {
//       newVectorData.push_back(jetToFbs(v));
//     }

//     auto fbsVectorData = fbs::CreateVectordPbdData3(
//         *builder, builder->CreateVectorOfStructs(newVectorData.data(), newVectorData.size()));
//     vectorDataList.push_back(fbsVectorData);
//   }
//   auto fbsVectorDataList = builder->CreateVector(vectorDataList);

//   // Copy the searcher
//   *fbsPbdData3 =
//       fbs::CreatePbdData3(*builder, _radius, _mass, _positionIdx, _velocityIdx, _forceIdx,
//                                 fbsScalarDataList, fbsVectorDataList);
// }

// void PbdData3::deserializedPbdSystemData(const fbs::PbdData3* fbsPbdData3) {
//   _scalarDataList.clear();
//   _vectorDataList.clear();

//   // Copy scalars
//   _radius = fbsPbdData3->radius();
//   _mass = fbsPbdData3->mass();
//   _positionIdx = static_cast<size_t>(fbsPbdData3->positionIdx());
//   _velocityIdx = static_cast<size_t>(fbsPbdData3->velocityIdx());
//   _forceIdx = static_cast<size_t>(fbsPbdData3->forceIdx());

//   // Copy data
//   auto fbsScalarDataList = fbsPbdData3->scalarDataList();
//   for (const auto& fbsScalarData : (*fbsScalarDataList)) {
//     auto data = fbsScalarData->data();

//     _scalarDataList.push_back(ScalarData(data->size()));

//     auto& newData = *(_scalarDataList.rbegin());

//     for (uint32_t i = 0; i < data->size(); ++i) {
//       newData[i] = data->Get(i);
//     }
//   }

//   auto fbsVectorDataList = fbsPbdData3->vectorDataList();
//   for (const auto& fbsVectorData : (*fbsVectorDataList)) {
//     auto data = fbsVectorData->data();

//     _vectorDataList.push_back(VectorData(data->size()));
//     auto& newData = *(_vectorDataList.rbegin());
//     for (uint32_t i = 0; i < data->size(); ++i) {
//       newData[i] = fbsToJet(*data->Get(i));
//     }
//   }

//   _numberOfParticles = _vectorDataList[0].size();
// }
