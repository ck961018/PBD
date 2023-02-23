#ifndef INCLUDE_JET_PBD_DATA3_H_
#define INCLUDE_JET_PBD_DATA3_H_

#include <jet/serialization.h>
#include <jet/vector2.h>
#include <jet/vector4.h>

#include <memory>
#include <vector>
#include <tuple>

#ifndef JET_DOXYGEN

namespace flatbuffers {

class FlatBufferBuilder;
template <typename T>
struct Offset;

}  // namespace flatbuffers

namespace jet {
namespace fbs {

struct PbdData3;

}
}  // namespace jet

#endif  // JET_DOXYGEN

namespace jet {

typedef std::tuple<size_t, size_t> IdxOfEdge;
typedef std::tuple<size_t, size_t, size_t, size_t> IdxOfTet;

class PbdData3 : public Serializable {
 public:
  //! Scalar data chunk.
  typedef Array1<double> ScalarData;

  //! Vector data chunk.
  typedef Array1<Vector3D> VectorData;

  //! Default constructor.
  PbdData3();

  //! Constructs particle system data with given number of particles.
  explicit PbdData3(size_t numberOfParticles);

  //! Copy constructor.
  PbdData3(const PbdData3& other);

  //! Destructor.
  virtual ~PbdData3();

  //!
  //! \brief      Resizes the number of particles of the container.
  //!
  //! 该处代码与原框架不同。PBD需要边与体的信息，部分信息存入了_scalarDataList，
  //! 由于边与体的数量可能大于点，原框架中resize同步更改_scalarDataList大小的方
  //! 案需要修改。
  //!
  //! \param[in]  newNumberOfParticles    New number of particles.
  //!
  void resize(size_t newNumberOfParticles);

  //! Returns the number of particles.
  size_t numberOfParticles() const;

  //!
  //! \brief      Adds a scalar data layer and returns its index.
  //!
  //! This function adds a new scalar data layer to the system. It can be used
  //! for adding a scalar attribute, such as temperature, to the particles.
  //!
  //! \params[in] initialVal  Initial value of the new scalar data.
  //!
  size_t addScalarData(double initialVal = 0.0);

  //!
  //! \brief      Adds a vector data layer and returns its index.
  //!
  //! This function adds a new vector data layer to the system. It can be used
  //! for adding a vector attribute, such as vortex, to the particles.
  //!
  //! \params[in] initialVal  Initial value of the new vector data.
  //!
  size_t addVectorData(const Vector3D& initialVal = Vector3D());

  //! Returns the mass of the particles.
  double mass() const;

  //! Sets the mass of the particles.
  virtual void setMass(double newMass);

  //! Returns the position array (immutable).
  ConstArrayAccessor1<Vector3D> positions() const;

  //! Returns the position array (mutable).
  ArrayAccessor1<Vector3D> positions();

  //! Returns the velocity array (immutable).
  ConstArrayAccessor1<Vector3D> velocities() const;

  //! Returns the velocity array (mutable).
  ArrayAccessor1<Vector3D> velocities();

  //! Returns the force array (immutable).
  ConstArrayAccessor1<Vector3D> forces() const;

  //! Returns the force array (mutable).
  ArrayAccessor1<Vector3D> forces();

  ConstArrayAccessor1<double> restLen() const;

  ArrayAccessor1<double> restLen();

  ConstArrayAccessor1<double> restVol() const;

  ArrayAccessor1<double> restVol();

  ConstArrayAccessor1<double> invMass() const;

  ArrayAccessor1<double> invMass();

  //! Returns custom scalar data layer at given index (immutable).
  ConstArrayAccessor1<double> scalarDataAt(size_t idx) const;

  //! Returns custom scalar data layer at given index (mutable).
  ArrayAccessor1<double> scalarDataAt(size_t idx);

  //! Returns custom vector data layer at given index (immutable).
  ConstArrayAccessor1<Vector3D> vectorDataAt(size_t idx) const;

  //! Returns custom vector data layer at given index (mutable).
  ArrayAccessor1<Vector3D> vectorDataAt(size_t idx);

  IdxOfEdge positionIdxOfEdge(size_t idx);

  ArrayAccessor1<IdxOfEdge> positionIdxOfEdges();

  ConstArrayAccessor1<IdxOfEdge> positionIdxOfEdges() const;

  IdxOfTet positionIdxOfTet(size_t idx);

  ArrayAccessor1<IdxOfTet> positionIdxOfTets();

  ConstArrayAccessor1<IdxOfTet> positionIdxOfTets() const;

  double lengthOfEdge(size_t idx);

  double volumeOfTet(size_t idx);

  //!
  //! \brief      Adds a particle to the data structure.
  //!
  //! This function will add a single particle to the data structure. For
  //! custom data layers, zeros will be assigned for new particles.
  //! However, this will invalidate neighbor searcher and neighbor lists. It
  //! is users responsibility to call
  //! ParticleSystemData3::buildNeighborSearcher and
  //! ParticleSystemData3::buildNeighborLists to refresh those data.
  //!
  //! \param[in]  newPosition The new position.
  //! \param[in]  newVelocity The new velocity.
  //! \param[in]  newForce    The new force.
  //!
  void addParticle(const Vector3D& newPosition, const Vector3D& newVelocity = Vector3D(),
                   const Vector3D& newForce = Vector3D());

  //!
  //! \brief      Adds particles to the data structure.
  //!
  //! This function will add particles to the data structure. For custom data
  //! layers, zeros will be assigned for new particles. However, this will
  //! invalidate neighbor searcher and neighbor lists. It is users
  //! responsibility to call ParticleSystemData3::buildNeighborSearcher and
  //! ParticleSystemData3::buildNeighborLists to refresh those data.
  //!
  //! \param[in]  newPositions  The new positions.
  //! \param[in]  newVelocities The new velocities.
  //! \param[in]  newForces     The new forces.
  //!
  void addParticles(
      const ConstArrayAccessor1<Vector3D>& newPositions,
      const ConstArrayAccessor1<Vector3D>& newVelocities = ConstArrayAccessor1<Vector3D>(),
      const ConstArrayAccessor1<Vector3D>& newForces = ConstArrayAccessor1<Vector3D>());

  void addEdge(const IdxOfEdge& positionIdxOfEdge);

  void addTet(const IdxOfTet& positionIdxOfTet);

  //   ! Serializes this particle system data to the buffer.
  void serialize(std::vector<uint8_t>* buffer) const override;

  //   ! Deserializes this particle system data from the buffer.
  void deserialize(const std::vector<uint8_t>& buffer) override;

  //! Copies from other particle system data.
  void set(const PbdData3& other);

  //! Copies from other particle system data.
  PbdData3& operator=(const PbdData3& other);

 protected:
  //   void serializeParticleSystemData(
  //       flatbuffers::FlatBufferBuilder* builder,
  //       flatbuffers::Offset<fbs::PbdData3>* fbsPbdSystemData) const;

  //   void deserializeParticleSystemData(const fbs::PbdData3* fbsPbdSystemData);

 private:
  // 单位体积质量
  double _mass = 1;
  size_t _numberOfParticles = 0;
  size_t _positionIdx;
  size_t _velocityIdx;
  size_t _forceIdx;

  // 这个mass是每个点上所有四面体体积*_mass/4的累加
  size_t _invMassIdx;
  size_t _restLenIdx;
  size_t _restVolIdx;

  std::vector<ScalarData> _scalarDataList;
  std::vector<VectorData> _vectorDataList;

  Array1<IdxOfTet> _positionIdxOfTets;
  Array1<IdxOfEdge> _positionIdxOfEdges;
};

//! Shared pointer type of PbdData3.
typedef std::shared_ptr<PbdData3> PbdData3Ptr;

}  // namespace jet

#endif  // INCLUDE_JET_PBD_DATA3_H_
