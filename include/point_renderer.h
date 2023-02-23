#ifndef POINT_RENDERER_H
#define POINT_RENDERER_H

#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <vector>

#include "texture.h"
#include "shader.h"

class PointRenderer {
 public:
  // Constructor (inits shaders/shapes)
  PointRenderer(Shader &shader);
  // Destructor
  ~PointRenderer();
  // Renders a defined quad textured with given sprite
  void draw(std::vector<glm::vec3> &position, float size = 1.0, glm::vec3 color = glm::vec3(1.0f));

 private:
  // Render state
  Shader shader;
  GLuint VAO, VBO;
  // Initializes and configures the quad's buffer and vertex attributes
  void initRenderData();
};

#endif