#include "Window.h"
#include <iostream>

static const float SCR_WIDTH_F = Core::SCR_WIDTH;
static const float SCR_HEIGHT_F = Core::SCR_HEIGHT;
static const char WINDOW_NAME[] = "fem";
static bool keys[1024];

namespace Core
{
    void KeyCallback(GLFWwindow* window, int key, int scancode, int action, int mode);

    GLFWwindow *CreateWindow()
    {
        return glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, WINDOW_NAME, nullptr, nullptr);
    }

    void TerminateWindow() {
        glfwTerminate();
    }

    void InitGL()
    {
        //Инициализация GLFW
        glfwInit();
        //Настройка GLFW
        //Задается минимальная требуемая версия OpenGL.
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);    //Мажорная
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);    //Минорная
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); //Установка профайла для которого создается контекст
        glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);                      //Выключение возможности изменения размера окна
    }

    bool InitWindowSetting(GLFWwindow *window)
    {
        glfwMakeContextCurrent(window);
        glfwSetFramebufferSizeCallback(window, FramebufferSizeCallback);
        glfwSetKeyCallback(window, KeyCallback);
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
        glewExperimental = GL_TRUE;

        if (glewInit() != GLEW_OK) {
            std::cout << "Failed to initialize GLEW" << std::endl;
            return false;
        }

        //glEnable(GL_DEPTH_TEST);
        //glEnable(GL_CULL_FACE);

        return true;
    }

    void KeyCallback(GLFWwindow* window, int key, int scancode, int action, int mode)
    {
        if (action == GLFW_PRESS)
            keys[key] = true;
        else if (action == GLFW_RELEASE)
            keys[key] = false;

        if (action == GLFW_PRESS && key == GLFW_KEY_ESCAPE)
        {
            glfwSetWindowShouldClose(window, GL_TRUE);
        }
    }

    void FramebufferSizeCallback(GLFWwindow *window, int width, int height)
    {
        glViewport(0, 0, width, height);
    }
}