#include <iostream>
#include "Window.h"
#include <Generator.h>
#include "Config.h"
#include <Builder.h>
#include <DMat.h>
#include <cassert>
#include "MSG.h"
#include "Shader.h"

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


    // build and compile our shader program
    // ------------------------------------
    Shader ourShader("3.3.shader.vs", "3.3.shader.fs"); // you can name your shader files however you like





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

    //for (int i = 1; i < q.size(); i++) std::cout << q[i] << " ";
    std::cout << std::endl;

    std::vector<float> gridVisual;

    double maxX = -1.0;
    double maxY = -1.0;
    double maxV = -1.0;
    double minV = 1000000000000;

    for (auto& fe : grid->finiteElements) {
        for (auto node : fe->nodes)
        {
            maxX = std::max(node->p.x, maxX);
        }
    }

    for (auto& fe : grid->finiteElements) {
        for (auto node : fe->nodes)
        {
            maxY = std::max(node->p.y, maxY);
        }
    }

    for (auto& fe : grid->finiteElements) {
        for (auto node : fe->nodes)
        {
            maxV = std::max(node->p.x, maxV);
        }
    }


    for (auto& fe : grid->finiteElements) {
        for (auto node : fe->nodes)
        {
            minV = std::min(node->p.y, minV);
        }
    }

    auto&& GetColor = [minV, maxV](double t)
    {
        return (t - minV) / maxV;
    };


    auto&& func = [&gridVisual, maxY, maxX, &GetColor](std::vector<Point> rectangle, std::vector<double> valuesOnRectangle, int i)
    {
        gridVisual.push_back(rectangle[i].x / maxX - 0.5f);
        gridVisual.push_back(rectangle[i].y / maxY - 0.5f);
        gridVisual.push_back(0.0);
        auto color = GetColor(valuesOnRectangle[i]);
        gridVisual.push_back(color);
        gridVisual.push_back(0.f); // GetColor
        gridVisual.push_back(1-color);
    };


    for (auto& fe : grid->finiteElements) {
        std::vector<Point> rectangle1 = { fe->nodes[0]->p, fe->nodes[1]->p, fe->nodes[3]->p, fe->nodes[4]->p };
        std::vector<double> valuesOnRectangle1 = { q[fe->nodes[0]->globalID], q[fe->nodes[1]->globalID], q[fe->nodes[3]->globalID], q[fe->nodes[4]->globalID] };

        func(rectangle1, valuesOnRectangle1, 0);
        func(rectangle1, valuesOnRectangle1, 1);
        func(rectangle1, valuesOnRectangle1, 2);
        func(rectangle1, valuesOnRectangle1, 1);
        func(rectangle1, valuesOnRectangle1, 2);
        func(rectangle1, valuesOnRectangle1, 3);

         std::vector<Point> rectangle2 = { fe->nodes[1]->p, fe->nodes[2]->p, fe->nodes[4]->p, fe->nodes[5]->p };
         std::vector<double> valuesOnRectangle2 = { q[fe->nodes[1]->globalID], q[fe->nodes[2]->globalID], q[fe->nodes[4]->globalID], q[fe->nodes[5]->globalID] };

         func(rectangle2, valuesOnRectangle2, 0);
         func(rectangle2, valuesOnRectangle2, 1);
         func(rectangle2, valuesOnRectangle2, 2);
         func(rectangle2, valuesOnRectangle2, 1);
         func(rectangle2, valuesOnRectangle2, 2);
         func(rectangle2, valuesOnRectangle2, 3);

         std::vector<Point> rectangle3 = { fe->nodes[3]->p, fe->nodes[4]->p, fe->nodes[6]->p, fe->nodes[7]->p };
         std::vector<double> valuesOnRectangle3 = { q[fe->nodes[3]->globalID], q[fe->nodes[4]->globalID], q[fe->nodes[6]->globalID], q[fe->nodes[7]->globalID] };

         func(rectangle3, valuesOnRectangle3, 0);
         func(rectangle3, valuesOnRectangle3, 1);
         func(rectangle3, valuesOnRectangle3, 2);
         func(rectangle3, valuesOnRectangle3, 1);
         func(rectangle3, valuesOnRectangle3, 2);
         func(rectangle3, valuesOnRectangle3, 3);

         std::vector<Point> rectangle4 = { fe->nodes[4]->p, fe->nodes[5]->p, fe->nodes[7]->p, fe->nodes[8]->p };
         std::vector<double> valuesOnRectangle4 = { q[fe->nodes[4]->globalID], q[fe->nodes[5]->globalID], q[fe->nodes[7]->globalID], q[fe->nodes[8]->globalID] };

         func(rectangle4, valuesOnRectangle4, 0);
         func(rectangle4, valuesOnRectangle4, 1);
         func(rectangle4, valuesOnRectangle4, 2);
         func(rectangle4, valuesOnRectangle4, 1);
         func(rectangle4, valuesOnRectangle4, 2);
         func(rectangle4, valuesOnRectangle4, 3);
    }
    


    // set up vertex data (and buffer(s)) and configure vertex attributes
    // ------------------------------------------------------------------
    /*float vertices[] = {
        // positions         // colors
         0.5f, -0.5f, 0.0f,  1.0f, 0.0f, 0.0f,  // bottom right
        -0.5f, -0.5f, 0.0f,  0.0f, 1.0f, 0.0f,  // bottom left
         0.0f,  0.5f, 0.0f,  0.0f, 0.0f, 1.0f   // top 
    };*/

    unsigned int VBO, VAO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    // bind the Vertex Array Object first, then bind and set vertex buffer(s), and then configure vertex attributes(s).
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * gridVisual.size(), gridVisual.data(), GL_STATIC_DRAW);

    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    // color attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    // You can unbind the VAO afterwards so other VAO calls won't accidentally modify this VAO, but this rarely happens. Modifying other
    // VAOs requires a call to glBindVertexArray anyways so we generally don't unbind VAOs (nor VBOs) when it's not directly necessary.
    // glBindVertexArray(0);

    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);


    // render loop
    // -----------
    while (!glfwWindowShouldClose(window))
    {
        // input
        // -----
        //processInput(window);

        // render
        // ------
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        // render the triangle
        ourShader.use();
        glBindVertexArray(VAO);
        glDrawArrays(GL_TRIANGLES, 0, 3 * gridVisual.size());

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // optional: de-allocate all resources once they've outlived their purpose:
    // ------------------------------------------------------------------------
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);

    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();
    //return 0;
    



    

}



//int _sz = 500;
//
//std::map<int, std::vector<std::pair<int, double>>> _l;
//std::vector<double> v, dot(_sz + 1, 0);
//
//for (int i = 0; i <= _sz; i++) {
//    v.push_back(rand() % 100);
//}
//
//for (int i = 1; i <= _sz; i++) {
//    for (int j = 1; j <= _sz; j++) {
//        _l[i].push_back({ j, i + j });
//        dot[i] += (i + j) * v[j];
//    }
//}
//
//DMat<double> A = DMat(_l, _sz);
//for (int i = 1; i <= _sz; i++) {
//    for (int j = 1; j <= _sz; j++) {
//        assert(A.get(i, j) == i + j);
//    }
//}
//
//
//std::vector<double> dot1 = A.product(v);
//for (int i = 1; i <= _sz; i++) {
//    assert(dot1[i] == dot[i]);
//}
