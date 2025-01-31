cmake_minimum_required (VERSION 3.13.4)
project (energy_efficient_touch-localization)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_FLAGS "-O3")

add_library(model STATIC model.cpp)
target_link_libraries(model /usr/local/lib/libyaml-cpp.a)

add_executable(main main.cpp)
target_link_libraries(main model)
target_link_libraries(main /usr/local/lib/libyaml-cpp.a)



