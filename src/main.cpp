#include <iostream>
#include "Window.h"
#include <Generator.h>
#include "Config.h"
#include <Builder.h>
#include <DMat.h>
#include "MSG.h"
#include "Line.h"
#include "Mesh.h"

int main() {
    Core::InitGL();
    auto window = Core::CreateWindow();
    if (!window)
    {
        std::cout << "Failed to create window" << std::endl;
        Core::TerminateWindow();
        return -1;
    }

    if (!Core::InitWindowSetting(window))
    {
        return -1;
    }

    Generator g;
    g.generateGrid(50, 50, 1, 1);

    Builder builder = Builder();
    Grid* grid = builder.createGrid();
    std::cout << "Generated grid" << std::endl;

    grid->calculateLocalMatrices();
    std::cout << "Calculated local matrices" << std::endl;

    DMat<double> A = grid->buildGlobalMatrix();
    grid->calculateLocalMatrices();
    std::cout << "Built global matrix" << std::endl;

    std::vector<double> b = grid->buildGlobalRP();
    std::cout << "Built right part vector" << std::endl;

    double __eps = 1e-15;
    double __maxiter = 1000;

    std::cout << "Solving SLE, maxiter = " << __maxiter << std::endl;
    std::vector<double> q = MSG::solve(A, b, __eps, __maxiter);
    std::cout << "Solved SLE." << std::endl;

    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    Mesh mesh(grid, q);

    float offsetSize = -0.55;
    glm::vec3 offsetAxis(offsetSize, offsetSize, 0.0);

    glm::vec3 axisXStart(0,0,0);
    glm::vec3 axisXEnd(1 + 0.05,0,0);
    glm::vec3 axisYStart(0,0,0);
    glm::vec3 axisYEnd(0,1 + 0.05,0);

    Line testLineX(axisXStart + offsetAxis, axisXEnd + offsetAxis);
    Line testLineY(axisYStart + offsetAxis, axisYEnd + offsetAxis);

    //testLine.setColor({1,1,1});
    while (!glfwWindowShouldClose(window))
    {
        // input
        // -----
        //processInput(window);

        // render
        // ------
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        mesh.draw();
        testLineX.draw();
        testLineY.draw();

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();

    return 0;
}
