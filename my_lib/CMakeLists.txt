cmake_minimum_required(VERSION 3.20)
project(my_lib)

add_subdirectory(my)

if (WITH_TESTS)
    enable_testing()
    add_subdirectory(test)
endif ()

