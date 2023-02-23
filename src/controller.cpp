#include "controller.h"
#include "resource_manager.h"
#include "point_renderer.h"
#include "triangle_renderer.h"
#include "skybox_renderer.h"

#include "pbd.hpp"
#include "camera.hpp"
#include "obj_constant.h"

#include <unordered_map>

PointRenderer* pointRenderer;
TriangleRenderer* triangleRenderer;
SkyboxRenderer* skyboxRenderer;

Camera camera(glm::vec3(0.0f, 0.0f, 3.0f));
std::ofstream logFile;

int numberOfFrames = 300000;
double fps = 60.0;
Frame frame(0, 1.0 / fps);

TriangleMesh3 mesh;

glm::vec3 lightPositions[] = {
    glm::vec3(-0.1f, 0.1f, 0.1f),
    glm::vec3(0.0f, 0.2f, 0.2f),
    glm::vec3(-0.5f, -0.2f, 0.1f),
    glm::vec3(0.1f, -0.3f, 0.4f),
};

unsigned int loadCubemap(std::vector<std::string> faces);

Controller::Controller(GLuint width, GLuint height)
    : State(CONTROLLER_ACTIVE), Keys(), Width(width), Height(height) {}

Controller::~Controller() {
  delete pointRenderer;
  delete triangleRenderer;
  delete skyboxRenderer;
}

void Controller::Init() {
  // 加载着色器
  ResourceManager::LoadShader("shaders/point.vs", "shaders/point.fs", nullptr, "point");
  ResourceManager::LoadShader("shaders/pbr.vs", "shaders/pbr.fs", nullptr, "triangle");
  ResourceManager::LoadShader("shaders/skybox.vs", "shaders/skybox.fs", nullptr, "skybox");

  // glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  glm::mat4 projection =
      glm::perspective(glm::radians(camera.Zoom), (float)Width / (float)Height, 0.1f, 100.0f);
  glm::mat4 view = camera.GetViewMatrix();
  glm::mat4 model = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, 0.0f));

  ResourceManager::GetShader("point").Use();
  ResourceManager::GetShader("point").SetMatrix4("projection", projection);
  ResourceManager::GetShader("point").SetMatrix4("view", view);
  ResourceManager::GetShader("point").SetMatrix4("model", model);

  ResourceManager::GetShader("triangle").Use();
  ResourceManager::GetShader("triangle").SetMatrix4("projection", projection);
  ResourceManager::GetShader("triangle").SetMatrix4("view", view);
  ResourceManager::GetShader("triangle").SetMatrix4("model", model);
  ResourceManager::GetShader("triangle").SetVector3f("camPos", camera.Position);
  ResourceManager::GetShader("triangle").SetInteger("skybox", 0);

  glm::vec3 lightColors[] = {glm::vec3(200.0f, 300.0f, 180.0f), glm::vec3(100.0f, 140.0f, 300.0f),
                             glm::vec3(100.0f, 200.0f, 210.0f), glm::vec3(300.0f, 300.0f, 300.0f)};
  for (size_t i = 0; i < 4; i++) {
    ResourceManager::GetShader("triangle")
        .SetVector3f(("lightColors[" + std::to_string(i) + "]").c_str(), lightColors[i]);
  }

  ResourceManager::GetShader("triangle").SetVector3f("albedo", {0.5f, 0.0f, 0.0f});
  ResourceManager::GetShader("triangle").SetFloat("ao", 1.0f);
  ResourceManager::GetShader("triangle").SetFloat("metallic", 1.0f);
  ResourceManager::GetShader("triangle").SetFloat("roughness", 1.0f);

  view = glm::mat4(glm::mat3(camera.GetViewMatrix()));
  ResourceManager::GetShader("skybox").Use();
  ResourceManager::GetShader("skybox").SetMatrix4("view", view);
  ResourceManager::GetShader("skybox").SetMatrix4("projection", projection);
  ResourceManager::GetShader("skybox").SetInteger("skybox", 0);

  std::vector<std::string> faces{
      R"(resources/textures/skybox/right.jpg)", R"(resources/textures/skybox/left.jpg)",
      R"(resources/textures/skybox/top.jpg)",   R"(resources/textures/skybox/bottom.jpg)",
      R"(resources/textures/skybox/front.jpg)", R"(resources/textures/skybox/back.jpg)"};
  GLuint cubemapTexture = loadCubemap(faces);

  auto numParticles = sizeof(Bunny::verts) / sizeof(Bunny::verts[0]) / 3;
  auto numTets = sizeof(Bunny::tetIds) / sizeof(Bunny::tetIds[0]) / 4;
  auto numEdges = sizeof(Bunny::tetEdgeIds) / sizeof(Bunny::tetEdgeIds[0]) / 2;
  auto numTriangles = sizeof(Bunny::tetSurfaceTriIds) / sizeof(Bunny::tetSurfaceTriIds[0]) / 3;

  for (int i = 0; i < numParticles * 3; i += 3) {
    mesh.addPoint({Bunny::verts[i], Bunny::verts[i + 1] + 1.0, Bunny::verts[i + 2]});
  }

  std::vector<std::tuple<size_t, size_t>> edges;
  std::vector<std::tuple<size_t, size_t, size_t, size_t>> tets;

  for (int i = 0; i < numEdges * 2; i += 2) {
    edges.push_back({Bunny::tetEdgeIds[i], Bunny::tetEdgeIds[i + 1]});
  }

  for (int i = 0; i < numTets * 4; i += 4) {
    tets.push_back(
        {Bunny::tetIds[i], Bunny::tetIds[i + 1], Bunny::tetIds[i + 2], Bunny::tetIds[i + 3]});
  }

  struct hash {
    size_t operator()(const Vector3D& e) const {
      return (static_cast<unsigned int>(e.x) * 73856093) ^
             (static_cast<unsigned int>(e.y) * 19349663) ^
             (static_cast<unsigned int>(e.z) * 83492791);
    }
  };
  std::unordered_map<Vector3D, size_t, hash> normalMap;
  for (int i = 0; i < numTriangles * 3; i += 3) {
    Point3UI pointIndex, normalIndex;
    pointIndex = {Bunny::tetSurfaceTriIds[i + 0], Bunny::tetSurfaceTriIds[i + 1],
                  Bunny::tetSurfaceTriIds[i + 2]};
    Vector3D A = mesh.point(pointIndex.x);
    Vector3D B = mesh.point(pointIndex.y);
    Vector3D C = mesh.point(pointIndex.z);
    auto normal = (A - B).cross(C - B).normalized();
    if (normalMap.find(normal) == normalMap.end()) {
      auto tmp = mesh.numberOfNormals();
      normalMap[normal] = tmp;
      mesh.addNormal(normal);

      normalIndex = {tmp, tmp, tmp};
    } else {
      auto tmp = normalMap[normal];

      normalIndex = {tmp, tmp, tmp};
    }
    mesh.addPointNormalTriangle(pointIndex, normalIndex);
  }

  pointRenderer = new PointRenderer(ResourceManager::GetShader("point"));
  triangleRenderer = new TriangleRenderer(ResourceManager::GetShader("triangle"),
                                          static_cast<GLsizei>(numTriangles), cubemapTexture);
  skyboxRenderer = new SkyboxRenderer(ResourceManager::GetShader("skybox"), cubemapTexture);

  // 设置log输出流
  std::string logFilename = R"(log.txt)";
  logFile.open(logFilename.c_str());
  if (logFile) {
    Logging::setAllStream(&logFile);
  }

  PBD::init(mesh, edges, tets);
}

void Controller::Update(GLfloat dt) {
  if (frame.index == numberOfFrames) {
    return;
  }
  PBD::update(mesh, frame);
  frame++;

  auto trans =
      glm::rotate(glm::mat4(1.0), glm::radians(20.0f * frame.index), glm::vec3(1.0, 0.5, 1.0));

  for (size_t i = 0; i < 4; i++) {
    ResourceManager::GetShader("triangle")
        .SetVector3f(("lightPositions[" + std::to_string(i) + "]").c_str(),
                     glm::vec4(lightPositions[i], 1.0) * trans);
  }
}

void Controller::ProcessInput(GLfloat dt) {
  if (this->Keys[GLFW_KEY_W]) {
    camera.move(UP, dt);
  }
  if (this->Keys[GLFW_KEY_S]) {
    camera.move(DOWN, dt);
  }
  if (this->Keys[GLFW_KEY_A]) {
    camera.move(LEFT, dt);
  }
  if (this->Keys[GLFW_KEY_D]) {
    camera.move(RIGHT, dt);
  }

  if (this->KeysMouse[GLFW_MOUSE_BUTTON_LEFT]) {
    auto dis = camera.Front;
    dis += glm::vec3{0, 0, 3};
    camera.Position += dis;
    auto offset = NowMousePos - OldMousePos;
    camera.rotate(-offset.x, offset.y);
    camera.Position -= dis;
  }
  if (this->KeysMouse[GLFW_MOUSE_BUTTON_RIGHT]) {
    auto offset = NowMousePos - OldMousePos;
    if (offset.x < 0.0) {
      camera.move(FORWARD, dt);
    } else {
      camera.move(BACKWARD, dt);
    }
  }

  glm::mat4 projection =
      glm::perspective(glm::radians(camera.Zoom), (float)Width / (float)Height, 0.1f, 100.0f);
  glm::mat4 view = camera.GetViewMatrix();
  glm::mat4 model = glm::mat4(1.0f);

  ResourceManager::GetShader("point").Use();
  ResourceManager::GetShader("point").SetMatrix4("projection", projection);
  ResourceManager::GetShader("point").SetMatrix4("view", view);
  ResourceManager::GetShader("point").SetMatrix4("model", model);

  ResourceManager::GetShader("triangle").Use();
  ResourceManager::GetShader("triangle").SetMatrix4("projection", projection);
  ResourceManager::GetShader("triangle").SetMatrix4("view", view);
  ResourceManager::GetShader("triangle").SetMatrix4("model", model);
  ResourceManager::GetShader("triangle").SetVector3f("camPos", camera.Position);

  view = glm::mat4(glm::mat3(camera.GetViewMatrix()));
  ResourceManager::GetShader("skybox").Use();
  ResourceManager::GetShader("skybox").SetMatrix4("view", view);
  ResourceManager::GetShader("skybox").SetMatrix4("projection", projection);
}

void Controller::Render() {
  std::vector<glm::vec3> particles;
  // tmp->positions().forEach([&](const auto& p) { particles.push_back({p.x, p.y, p.z}); });
  // pointRenderer->draw(particles, 3.0, {0.0, 1.0, 0.0});

  particles.clear();
  for (int i = 0; i < mesh.numberOfPoints(); i++) {
    auto _point = mesh.point(i);
    particles.push_back({_point.x, _point.y, _point.z});
  }

  std::vector<glm::vec3> normal;
  for (int i = 0; i < mesh.numberOfNormals(); i++) {
    auto _normal = mesh.normal(i);
    normal.push_back({_normal.x, _normal.y, _normal.z});
  }

  std::vector<glm::ivec3> positionIndex, normalIndex;
  for (int i = 0; i < mesh.numberOfTriangles(); i++) {
    auto _positionIndex = mesh.pointIndex(i);
    auto _normalIndex = mesh.normalIndex(i);
    positionIndex.push_back({_positionIndex.x, _positionIndex.y, _positionIndex.z});
    normalIndex.push_back({_normalIndex.x, _normalIndex.y, _normalIndex.z});
  }

  // printf("tris: %d\n", mesh.numberOfTriangles());

  triangleRenderer->draw(particles, positionIndex, normal, normalIndex);

  skyboxRenderer->draw();
}

unsigned int loadCubemap(std::vector<std::string> faces) {
  unsigned int textureID;
  glGenTextures(1, &textureID);
  glBindTexture(GL_TEXTURE_CUBE_MAP, textureID);

  int width, height, nrChannels;
  for (unsigned int i = 0; i < faces.size(); i++) {
    // printf("%s\n", faces[i].c_str());
    unsigned char* data = stbi_load(faces[i].c_str(), &width, &height, &nrChannels, 0);
    if (data) {
      glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, GL_RGB, width, height, 0, GL_RGB,
                   GL_UNSIGNED_BYTE, data);
      stbi_image_free(data);
    } else {
      printf("wrong image: %I32d", i);
      stbi_image_free(data);
    }
  }
  glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

  return textureID;
}