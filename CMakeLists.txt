cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(solver_fluidos VERSION 1.0)

# Flags de compilación adicionales (útil para depurar)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -g")

add_executable(main config_malla.hpp config_CF.hpp malla_por_bloques.cpp malla_por_bloques.hpp
        condiciones_de_frontera.cpp
        condiciones_de_frontera.hpp
        main.cpp
        esquemas_de_discretizacion.hpp
        config_control.hpp
        solvers_lineales.cpp
        solvers_lineales.hpp
        calculo_del_error.cpp
        calculo_del_error.hpp
        escritura.cpp
        escritura.hpp
        solvers_lineales.cpp
)
