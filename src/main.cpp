#include <iostream>
#include "Window.h"
#include <Generator.h>
#include "Config.h"
#include <Builder.h>
#include <DMat.h>
#include <cassert>
#include "MSG.h"


int main() {
    //Core::InitGL();
    //auto window = Core::CreateWindow();
    //if (!window)
    //{
    //    std::cout << "Failed to create window" << std::endl;
    //    Core::TerminateWindow();
    //    return -1;
    //}

    //if (!Core::InitWindowSetting(window))
    //{
    //    return -1;
    //}

    //while (!glfwWindowShouldClose(window))
    //{
    //    glfwPollEvents();
    //    glClearColor(0.45f, 0.58f, 0.68f, 1.0f);
    //    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    //    // Paste code for visual here

    //    glfwSwapBuffers(window);
    //}

    //glfwTerminate();




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

    std::cout << "Solving SLE, maxiter = " << __maxiter  << std::endl;
    std::vector<double> q = MSG::solve(A, b, __eps, __maxiter);
    std::cout << "Solved SLE." << std::endl;

    //for (int i = 1; i < q.size(); i++) std::cout << q[i] << " ";
    std::cout << std::endl;



    delete grid;

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
