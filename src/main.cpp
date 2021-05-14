#include <iostream>
#include "Window.h"

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

    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();
        glClearColor(0.45f, 0.58f, 0.68f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Paste code for visual here

        glfwSwapBuffers(window);
    }

    glfwTerminate();
    return 0;
}
