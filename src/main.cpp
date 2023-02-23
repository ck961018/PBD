#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "controller.h"
#include "resource_manager.h"

#include <iostream>

// GLFW function declarations
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode);
void mouse_callback(GLFWwindow* window, double xposIn, double yposIn);
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);

// The Width of the screen
const unsigned int SCREEN_WIDTH = 800;
// The height of the screen
const unsigned int SCREEN_HEIGHT = 600;

Controller Sph3D(SCREEN_WIDTH, SCREEN_HEIGHT);
bool firstMouse = true;

int main(int argc, char* argv[]) {
  glfwInit();
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif
  glfwWindowHint(GLFW_RESIZABLE, false);

  GLFWwindow* window = glfwCreateWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "Sph3D", nullptr, nullptr);
  glfwMakeContextCurrent(window);

  // glad: load all OpenGL function pointers
  // ---------------------------------------
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }

  glfwSetKeyCallback(window, key_callback);
  glfwSetCursorPosCallback(window, mouse_callback);
  glfwSetMouseButtonCallback(window, mouse_button_callback);
  glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

  // OpenGL configuration
  // --------------------
  glViewport(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT);

  glEnable(GL_DEPTH_TEST);

  // initialize game
  // ---------------
  Sph3D.Init();

  // deltaTime variables
  // -------------------
  float deltaTime = 0.0f;
  float lastFrame = 0.0f;

  Sph3D.State = CONTROLLER_ACTIVE;

  glfwSwapInterval(1);

  while (!glfwWindowShouldClose(window)) {
    // calculate delta time
    // --------------------
    float currentFrame = static_cast<float>(glfwGetTime());
    deltaTime = currentFrame - lastFrame;
    lastFrame = currentFrame;
    glfwPollEvents();

    // manage user input
    // -----------------
    Sph3D.ProcessInput(deltaTime);

    Sph3D.Update(deltaTime);

    // render
    // ------
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    Sph3D.Render();

    glfwSwapBuffers(window);
  }

  // delete all resources as loaded using the resource manager
  // ---------------------------------------------------------
  ResourceManager::Clear();

  glfwTerminate();
  return 0;
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode) {
  // when a user presses the escape key, we set the WindowShouldClose property to true, closing the
  // application
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) glfwSetWindowShouldClose(window, true);
  if (key >= 0 && key < 1024) {
    if (action == GLFW_PRESS)
      Sph3D.Keys[key] = true;
    else if (action == GLFW_RELEASE) {
      Sph3D.Keys[key] = false;
    }
  }
}

void mouse_callback(GLFWwindow* window, double xposIn, double yposIn) {
  float xpos = static_cast<float>(xposIn);
  float ypos = static_cast<float>(yposIn);

  if (firstMouse) {
    Sph3D.OldMousePos = {xpos, ypos};
    firstMouse = false;
  } else {
    Sph3D.OldMousePos = Sph3D.NowMousePos;
  }
  Sph3D.NowMousePos = {xpos, ypos};
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
  if (action == GLFW_PRESS) {
    Sph3D.KeysMouse[button] = true;
  } else if (action == GLFW_RELEASE) {
    Sph3D.KeysMouse[button] = false;
  }
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
  // make sure the viewport matches the new window dimensions; note that width and
  // height will be significantly larger than specified on retina displays.
  glViewport(0, 0, width, height);
}