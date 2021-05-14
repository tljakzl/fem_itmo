cmake_minimum_required(VERSION 3.19)
project(fem_itmo)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(LIBSDIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(SRCDIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(BINDIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# init exe
file(GLOB fem_SRC
        "${SRCDIR}/*.h"
        "${SRCDIR}/*.cpp"
        )

add_executable(fem_itmo ${fem_SRC} src/Window.cpp src/Window.h)

# set output directory exe
if ( MSVC )
    set_target_properties(fem_itmo PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG ${BINDIR})
    set_target_properties(fem_itmo PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE ${BINDIR})
    set_target_properties(fem_itmo PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${BINDIR})
else()
    set_target_properties(fem_itmo PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG ${BINDIR})
endif ( MSVC )

if( MSVC )
    set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT "fem_itmo")
    set_target_properties(fem_itmo PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${BINDIR}")
endif( MSVC )

# add OpenGl
find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIRS})
if(WIN32)
    target_compile_options(fem_itmo PRIVATE /EHsc)
    target_link_libraries(fem_itmo ${OPENGL_LIBRARY})
else()
    target_link_libraries(fem_itmo GL)
endif()

# add source
target_include_directories(fem_itmo PUBLIC ${SRCDIR})

# add glfw
target_include_directories(fem_itmo PRIVATE ${LIBSDIR}/glfw/include)
add_subdirectory(${LIBSDIR}/glfw)
target_link_libraries(fem_itmo glfw)

# add GLEW
target_include_directories(fem_itmo PRIVATE ${LIBSDIR}/glew-cmake/include)
add_subdirectory(${LIBSDIR}/glew-cmake)
target_link_libraries(fem_itmo libglew_static)

# add glm
target_include_directories(fem_itmo PRIVATE ${LIBSDIR}/glm)
add_subdirectory(${LIBSDIR}/glm)