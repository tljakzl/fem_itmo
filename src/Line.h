#pragma once
#include <vector>
#include <glm/glm.hpp>
#include "Shader.h"

class Line {
public:
    explicit Line(glm::vec3 start, glm::vec3 end);
    void Init();
    ~Line();
    void draw();

private:
    Line() = delete;
public:
    unsigned int VBO;
    unsigned int VAO;
    glm::vec3 startPoint;
    glm::vec3 endPoint;
    glm::vec3 lineColor;
    std::vector<float> vertices;
    Shader shader;
};