cmake_minimum_required(VERSION 3.12.0)
project(Ptero LANGUAGES CXX C)

include(../common/cmake-common.cmake)

add_executable(navier_semi_implicit src/main.cpp src/Navier2.cpp)
deal_ii_setup_target(navier_semi_implicit)