#include "Line.h"

/*int setMVP(mat4 mvp) {
    MVP = mvp;
}

int setColor(glm::vec3 color) {
    lineColor = color;
}*/

Line::~Line() {
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
}

void Line::draw() {
    shader.use();
    shader.setVec3("color", lineColor);
    //glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "MVP"), 1, GL_FALSE, &MVP[0][0]);

    glBindVertexArray(VAO);
    glDrawArrays(GL_LINES, 0, vertices.size()/3);
    glBindVertexArray(0);
}

Line::Line(glm::vec3 start, glm::vec3 end)
: VBO(0)
, VAO(0)
, startPoint(start)
, endPoint(end)
, lineColor(glm::vec3(1, 1, 1))
, vertices({start.x, start.y, start.z, end.x, end.y, end.z,})
, shader("lineShader.vs", "lineShader.fs")
{
    Init();
}

void Line::Init() {
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * vertices.size(), vertices.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void *) 0);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}


