cmake_minimum_required(VERSION 3.18)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(gaussian_duct)

# find sources
file(GLOB SRCFILES *.cpp)
file(GLOB HPPFILES *.hpp)

add_executable(${PROJECT_NAME} ${SRCFILES} ${HPPFILES})