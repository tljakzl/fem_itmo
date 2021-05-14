#pragma once

#include <GL/glew.h>
#include <GLFW/glfw3.h>

namespace Core
{
    static const unsigned int SCR_WIDTH  = 16 * 100;
    static const unsigned int SCR_HEIGHT = 9 * 100;

    GLFWwindow* CreateWindow();
    void InitGL();
    void TerminateWindow();
    bool InitWindowSetting(GLFWwindow* window);
    void FramebufferSizeCallback(GLFWwindow* window, int width, int height);
}



