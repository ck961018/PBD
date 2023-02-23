#include "triangle_renderer.h"

TriangleRenderer::TriangleRenderer(Shader &shader, GLsizei trianglesNum, GLuint _cubemapTexture) {
  this->shader = shader;
  this->cubemapTexture = _cubemapTexture;
  this->initRenderData(trianglesNum);
}

TriangleRenderer::~TriangleRenderer() {
  glDeleteVertexArrays(1, &VAO);
  glDeleteBuffers(1, &VBO);
  glDeleteBuffers(1, &EBO);
}

void TriangleRenderer::draw(std::vector<glm::vec3> &position,
                            std::vector<glm::ivec3> &positionIndex, std::vector<glm::vec3> &normal,
                            std::vector<glm::ivec3> &normalIndex) {
  // prepare transformations
  this->shader.Use();
  this->shader.SetVector3f("spriteColor", {0.2, 0.5, 0.3});

  glBindVertexArray(VAO);
  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  glBufferSubData(GL_ARRAY_BUFFER, 0, position.size() * sizeof(glm::vec3), position.data());
  glBufferSubData(GL_ARRAY_BUFFER, position.size() * sizeof(glm::vec3),
                  normal.size() * sizeof(glm::vec3), normal.data());

  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 1 * sizeof(glm::vec3), nullptr);
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 1 * sizeof(glm::vec3),
                        (void *)(position.size() * sizeof(glm::vec3)));

  // glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(position.size()));

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
  glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, positionIndex.size() * sizeof(glm::ivec3),
                  positionIndex.data());

  glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, positionIndex.size() * sizeof(glm::ivec3),
                  normalIndex.size() * sizeof(glm::ivec3), normalIndex.data());

  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_CUBE_MAP, cubemapTexture);
  glDrawElements(GL_TRIANGLES, static_cast<GLsizei>(positionIndex.size() * 3), GL_UNSIGNED_INT, 0);

  glBindVertexArray(0);
}

void TriangleRenderer::initRenderData(GLsizei trianglesNum) {
  glGenVertexArrays(1, &VAO);
  glGenBuffers(1, &VBO);
  glGenBuffers(1, &EBO);

  glBindVertexArray(VAO);
  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  glBufferData(GL_ARRAY_BUFFER, trianglesNum * 3 * sizeof(glm::vec3), nullptr, GL_DYNAMIC_DRAW);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, trianglesNum * 3 * sizeof(glm::ivec3), nullptr,
               GL_DYNAMIC_DRAW);
}