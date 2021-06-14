#include <vector>
#include <glm/glm.hpp>
#include "Shader.h"
#include "Point.h"

class Grid;

class Mesh {
public:
    Mesh(Grid* grid, const std::vector<double>& q);
    void draw();
    ~Mesh();

private:
    void AddPoint(const Point& point, double value);
    void AddQuad(const std::vector<Point>& rectangle, const std::vector<double>& valuesOnRectangle);
    double NormalizeValue(double value, double min, double max);

private:
    unsigned int VBO;
    unsigned int VAO;
    Shader shader;
    std::vector<float> gridVisual;
    double maxV;
    double minV;
    double minY;
    double maxY;
    double minX;
    double maxX;
};

