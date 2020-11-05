cmake_minimum_required(VERSION 3.2.0)
project(quadprog_test)

SET(SRC_PATH 
    ${CMAKE_CURRENT_LIST_DIR}/src
    
) 

SET(EXAMPLES_PATH
    ${CMAKE_CURRENT_LIST_DIR}/examples
)

include_directories(
    ${SRC_PATH}
    ${EXAMPLES_PATH}
)

file(GLOB SRC_LIST 
    ${SRC_PATH}/*.c
    ${EXAMPLES_PATH}/*.c
)

# add_library(test_lib STATIC ${SRC_LIST})
add_executable(demo ${SRC_LIST})
target_link_libraries(demo m)