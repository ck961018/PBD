#ifndef CONTROLL_H
#define CONTROLL_H

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

enum ControllerState { CONTROLLER_ACTIVE, CONTROLLER_MENU, CONTROLLER_WIN };

class Controller {
 public:
  ControllerState State;
  GLboolean Keys[1024];
  GLboolean KeysMouse[12];
  glm::vec2 OldMousePos;
  glm::vec2 NowMousePos;

  GLuint Width, Height;
  Controller(GLuint width, GLuint height);
  ~Controller();
  void Init();
  void ProcessInput(GLfloat dt);
  void Update(GLfloat dt);
  void Render();
};

#endif