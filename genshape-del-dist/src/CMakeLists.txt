file(GLOB SRC_FILES *.cpp)
add_executable(${PROJECT_NAME} ${SRC_FILES})
# target_link_libraries(${PROJECT_NAME} PUBLIC igl::glfw )
# target_link_libraries(${PROJECT_NAME} PRIVATE igl_copyleft::tetgen )


