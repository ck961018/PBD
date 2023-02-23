#ifndef TRIANGLE_RENDERER_H
#define TRIANGLE_RENDERER_H

#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <vector>

#include "texture.h"
#include "shader.h"

class TriangleRenderer {
 public:
  // Constructor (inits shaders/shapes)
  TriangleRenderer(Shader &shader, GLsizei trianglesNum, GLuint cubemapTexture);
  // Destructor
  ~TriangleRenderer();
  // Renders a defined quad textured with given sprite
  void draw(std::vector<glm::vec3> &position, std::vector<glm::ivec3> &positionIndex,
            std::vector<glm::vec3> &normal, std::vector<glm::ivec3> &normalIndex);

 private:
  // Render state
  Shader shader;
  GLuint VAO, VBO, EBO, cubemapTexture;
  // Initializes and configures the quad's buffer and vertex attributes
  void initRenderData(GLsizei trianglesNum);
};

#endif