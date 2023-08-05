cmake_minimum_required(VERSION 3.12.0)
project(Ptero LANGUAGES CXX C)

include(./common/cmake-common.cmake)

add_executable(navier src/main.cpp src/Navier.cpp)
deal_ii_setup_target(navier)