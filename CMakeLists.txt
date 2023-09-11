cmake_minimum_required(VERSION 3.0.0)
project(
    efp-control
    VERSION 0.1
    DESCRIPTION "Eager functional programming library for C++ 11"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED On)
set(CMAKE_CXX_COMPILER /usr/bin/clang++)

# set(CMAKE_CXX_COMPILER /usr/bin/g++)
find_package(Catch2 3 REQUIRED)

include(CTest)
include(Catch)
enable_testing()

add_library(efp-control
    INTERFACE)

target_link_libraries(efp-control
    INTERFACE
    efp
    eigen)

target_include_directories(efp-control
    INTERFACE
    "./include")

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

add_subdirectory(lib)
add_subdirectory(test)
add_subdirectory(include)