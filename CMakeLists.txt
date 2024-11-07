cmake_minimum_required(VERSION 3.5.0)
project(efp_control VERSION 0.1.0 LANGUAGES CXX)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

# temp C++ 14 for auto lambda
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED true)

include("cmake/CPM.cmake")

CPMAddPackage(
    NAME efp
    VERSION 0.1.0-rc.6
    URL "https://github.com/cwahn/efp/archive/refs/tags/v0.1.0-rc.6.tar.gz"
)

CPMAddPackage(
    NAME eigen
    VERSION 3.4.0
    URL "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz"
)

file(GLOB EFP_CONTROL_SOURCE src/*.cpp)
add_library(efp_control STATIC ${EFP_CONTROL_SOURCE})

target_include_directories(efp_control
    PUBLIC
    "./include")

target_link_libraries(efp_control
    PUBLIC
    efp
    eigen)

# Add sub-directories only if this is the top-level project
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    add_subdirectory(test)
    add_subdirectory(example)
endif()