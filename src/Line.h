#pragma once
#include <vector>
#include <glm/glm.hpp>

class Line {
public:
    Line(glm::vec3 start, glm::vec3 end);
    int draw();
    ~Line();

private:
    int shaderProgram;
    unsigned int VBO, VAO;
    std::vector<float> vertices;
    glm::vec3 startPoint;
    glm::vec3 endPoint;
    //mat4 MVP;
    glm::vec3 lineColor;
};