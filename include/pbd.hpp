#include <jet/jet.h>

#include <tuple>

using namespace jet;

namespace PBD {

PbdSolver3Ptr solver;

void init(TriangleMesh3 &mesh, std::vector<std::tuple<size_t, size_t>> edges,
          std::vector<std::tuple<size_t, size_t, size_t, size_t>> tets) {
  solver = PbdSolver3::builder().makeShared();
  auto pbdData = solver->pbdData();
  pbdData->addParticles(mesh.points());
  // 加约束
  for (const auto &e : edges) {
    pbdData->addEdge(e);
  }

  for (const auto &t : tets) {
    pbdData->addTet(t);
  }

  auto plane3 = Plane3::builder().withNormal({0, 1, 0}).withPoint({0, -1, 0}).makeShared();
  auto collider = RigidBodyCollider3::builder().withSurface(plane3).makeShared();
  solver->setCollider(collider);

  solver->setNumberOfFixedSubTimeSteps(10);
  solver->setGravity({0, -1, 0});
}

void update(TriangleMesh3 &mesh, Frame &frame) {
  solver->update(frame);
  mesh.points().clear();
  auto positions = solver->pbdData()->positions();
  for (int i = 0; i < positions.size(); i++) {
    mesh.addPoint(positions[i]);
  }
}
}  // namespace PBD