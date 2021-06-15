#include "Mesh.h"
#include "Color.h"
#include "Grid.h"
#include <numeric>

double Mesh::NormalizeValue(double value, double min, double max) {
    return (value - min) / (max - min);
}

void Mesh::AddPoint(const Point& point, double value) {
    auto offset = -0.5f;
    auto x = NormalizeValue(point.x, minX, maxX);
    auto y = NormalizeValue(point.y, minY, maxY);
    auto z = 0.f;

    auto normValue = NormalizeValue(value, minV, maxV);
    auto color = GetColor(normValue);

    gridVisual.push_back(x + offset);
    gridVisual.push_back(y + offset);
    gridVisual.push_back(z);
    gridVisual.push_back(color.r());
    gridVisual.push_back(color.g());
    gridVisual.push_back(color.b());
}

void Mesh::AddQuad(const std::vector<Point>& rectangle, const std::vector<double>& valuesOnRectangle) {
    AddPoint(rectangle[0], valuesOnRectangle[0]);
    AddPoint(rectangle[1], valuesOnRectangle[1]);
    AddPoint(rectangle[2], valuesOnRectangle[2]);
    AddPoint(rectangle[1], valuesOnRectangle[1]);
    AddPoint(rectangle[2], valuesOnRectangle[2]);
    AddPoint(rectangle[3], valuesOnRectangle[3]);
}

Mesh::Mesh(Grid* grid, const std::vector<double>& q)
: VBO(0)
, VAO(0)
, shader("3.3.shader.vs", "3.3.shader.fs")
, gridVisual()
, maxV(std::numeric_limits<double>::min())
, minV(std::numeric_limits<double>::max())
, minY(std::numeric_limits<double>::max())
, maxY(std::numeric_limits<double>::min())
, minX(std::numeric_limits<double>::max())
, maxX(std::numeric_limits<double>::min())
{
    for (auto &fe : grid->finiteElements) {
        for (auto node : fe->nodes) {
            maxX = std::max(node->p.x, maxX);
            maxY = std::max(node->p.y, maxY);
            minX = std::min(node->p.x, minX);
            minY = std::min(node->p.y, minY);
            maxV = std::max(q[node->globalID], maxV);
            minV = std::min(q[node->globalID], minV);
        }
    }

    for (auto &fe : grid->finiteElements) {
        std::vector<Point> rectangle1 = {fe->nodes[0]->p, fe->nodes[4]->p, fe->nodes[8]->p, fe->nodes[12]->p};
        std::vector<double> valuesOnRectangle1 = {q[fe->nodes[0]->globalID] + q[fe->nodes[1]->globalID] + q[fe->nodes[2]->globalID] + q[fe->nodes[3]->globalID],
                                                  q[fe->nodes[4]->globalID] + q[fe->nodes[5]->globalID] + q[fe->nodes[6]->globalID] + q[fe->nodes[7]->globalID],
                                                  q[fe->nodes[8]->globalID] + q[fe->nodes[9]->globalID] + q[fe->nodes[10]->globalID] + q[fe->nodes[11]->globalID],
                                                  q[fe->nodes[12]->globalID] + q[fe->nodes[13]->globalID] + q[fe->nodes[14]->globalID] + q[fe->nodes[15]->globalID]};
        AddQuad(rectangle1, valuesOnRectangle1);
    }

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    // bind the Vertex Array Object first, then bind and set vertex buffer(s), and then configure vertex attributes(s).
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * gridVisual.size(), gridVisual.data(), GL_STATIC_DRAW);

    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *) 0);
    glEnableVertexAttribArray(0);
    // color attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void *) (3 * sizeof(float)));
    glEnableVertexAttribArray(1);
}

Mesh::~Mesh() {
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
}

void Mesh::draw() {
    shader.use();
    glBindVertexArray(VAO);
    glDrawArrays(GL_TRIANGLES, 0, gridVisual.size()/6);
    glBindVertexArray(0);
}
