#include "point_renderer.h"

PointRenderer::PointRenderer(Shader &shader) {
  this->shader = shader;
  this->initRenderData();
}

PointRenderer::~PointRenderer() {
  glDeleteVertexArrays(1, &VAO);
  glDeleteBuffers(1, &VBO);
}

void PointRenderer::draw(std::vector<glm::vec3> &position, float size, glm::vec3 color) {
  // prepare transformations
  this->shader.Use();

  this->shader.SetVector3f("spriteColor", color);

  glBindVertexArray(VAO);
  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  glBufferSubData(GL_ARRAY_BUFFER, 0, position.size() * sizeof(glm::vec3), position.data());
  glPointSize(size);  // 设置顶点大小

  glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(position.size()));
  glBindVertexArray(0);
}

void PointRenderer::initRenderData() {
  glGenVertexArrays(1, &VAO);
  glGenBuffers(1, &VBO);

  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  glBufferData(GL_ARRAY_BUFFER, 50000, nullptr, GL_STATIC_DRAW);

  glBindVertexArray(VAO);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void *)0);

  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
}