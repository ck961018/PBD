#ifndef SKYBOX_RENDERER_H
#define SKYBOX_RENDERER_H

#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <vector>

#include "texture.h"
#include "shader.h"
#include <stb_image.h>

class SkyboxRenderer {
 public:
  // Constructor (inits shaders/shapes)
  SkyboxRenderer(Shader &shader, GLuint _cubemapTexture);
  // Destructor
  ~SkyboxRenderer();
  // Renders a defined quad textured with given sprite
  void draw();

 private:
  // Render state
  Shader shader;
  GLuint VAO, VBO, cubemapTexture;
  // Initializes and configures the quad's buffer and vertex attributes
  void initRenderData();
};

#endif