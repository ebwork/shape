cmake_minimum_required(VERSION 3.16)
project(genshape-dd)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set( LIBIGL_COPYLEFT_TETGEN ON )
option(LIBIGL_COPYLEFT_TETGEN  "Build target igl_copyleft::tetgen"  ${LIBIGL_TOPLEVEL_PROJECT})

# Libigl
include(libigl)

# Enable the target igl::glfw
igl_include(glfw)

# Add your project files
add_subdirectory( src )
target_link_libraries(${PROJECT_NAME} PUBLIC igl::glfw )

# target_link_libraries(${PROJECT_NAME} PUBLIC igl::core igl::opengl_glfw)
# target_link_libraries(${PROJECT_NAME} PUBLIC igl::core igl::glfw)
# target_link_libraries(${PROJECT_NAME} PRIVATE igl_copyleft::tetgen )
# target_link_libraries(${PROJECT_NAME} PUBLIC igl::glfw igl::tetgen )
# target_link_libraries(${PROJECT_NAME} PUBLIC igl_copyleft::tetgen )
# target_link_libraries(${PROJECT_NAME} PUBLIC igl::tetgen )
if ( TARGET igl_copyleft::tetgen)
    target_link_libraries(${PROJECT_NAME} PUBLIC igl_copyleft::tetgen )
endif()

# ----------------------------------------


