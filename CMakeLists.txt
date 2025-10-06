cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(solver_fluidos VERSION 1.0)

# Flags de compilación adicionales (útil para depurar)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++23 -Wall -Wextra -g")

file(GLOB_RECURSE SRC_FILES src/*.cpp)
set(CMAKE_BUILD_TYPE Debug)
add_executable(main main.cpp ${SRC_FILES}
        src/correccion_campos.cpp
        include/correccion_campos.hpp)
target_include_directories(main PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
