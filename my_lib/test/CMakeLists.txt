cmake_minimum_required(VERSION 3.20)

enable_testing()
find_package(GTest REQUIRED)

file(GLOB_RECURSE TEST_SOURCES test_*.cpp)

foreach (TEST_FILE IN LISTS TEST_SOURCES)
    get_filename_component(TEST_FILE_NAME ${TEST_FILE} NAME_WE)
    set(EXECUTABLE_NAME run_${TEST_FILE_NAME})

    add_executable(${EXECUTABLE_NAME} ${TEST_FILE})
    
    add_subdirectory(../ dim/)
    target_link_libraries(${EXECUTABLE_NAME} PRIVATE my_lib GTest::gtest_main)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR
       CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        target_compile_options(${EXECUTABLE_NAME} PRIVATE -Werror -Wall -Wextra -pedantic -Wshadow -Wconversion -Wsign-conversion)
    endif()
    
    add_test(NAME ${EXECUTABLE_NAME} COMMAND ${EXECUTABLE_NAME})
endforeach ()
