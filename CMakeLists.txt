# Como correr fuera de CLion:

# MODO Debug:
# cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
# cmake --build build -j

# MODO Release:
# cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
# cmake --build build -j

# Modo RelWithDebInfo (con PARALLEL, -O2 y símbolos):
# cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
# cmake --build build -j

cmake_minimum_required(VERSION 3.10)
project(solver_fluidos VERSION 1.0 LANGUAGES CXX)

# Estándar C++23 sin extensiones GNU (usa -std=c++23 en lugar de -std=gnu++23)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# compile_commands.json para LSP
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Tipos de build: si no especificas, por defecto Debug
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "No build type selected, defaulting to Debug")
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

# Fuentes
file(GLOB_RECURSE SRC_FILES src/*.cpp)
add_executable(main main.cpp ${SRC_FILES}
        src/convergencia.cpp
        include/convergencia.hpp)
target_include_directories(main PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# ===== oneTBB (necesario para std::execution::par en libstdc++) =====
# Intento 1: paquete CONFIG moderno (TBB::tbb)
find_package(TBB QUIET CONFIG)
if(NOT TBB_FOUND)
  # Intento 2: modo clásico
  find_package(TBB REQUIRED)
endif()

# Determina el target correcto (TBB::tbb moderno o 'tbb' clásico)
set(_TBB_TARGET "")
if(TARGET TBB::tbb)
  set(_TBB_TARGET TBB::tbb)
elseif(TARGET tbb)
  set(_TBB_TARGET tbb)
endif()

if(_TBB_TARGET STREQUAL "")
  message(FATAL_ERROR "No se encontró un target de TBB exportado. Revisa que 'libtbb-dev' esté instalado.")
endif()

target_link_libraries(main PRIVATE ${_TBB_TARGET})

# ===== Definiciones y flags por configuración =====
# Define PARALLEL en Release/RelWithDebInfo; DEBUG en Debug
target_compile_definitions(main PRIVATE
  $<$<CONFIG:Debug>:DEBUG>
  $<$<CONFIG:Release>:NDEBUG;PARALLEL>
  $<$<CONFIG:RelWithDebInfo>:NDEBUG;PARALLEL>
)

# Warnings y sanitizers solo en Debug (GCC/Clang)
target_compile_options(main PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall;-Wextra;-Wpedantic>
  $<$<AND:$<CONFIG:Debug>,$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>>:-fsanitize=address,undefined>
)
target_link_options(main PRIVATE
  $<$<AND:$<CONFIG:Debug>,$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>>:-fsanitize=address,undefined>
)

# Optimizaciones por build type (puedes ajustar a tu gusto)
# Nota: CMake ya añade -O3/-O2 y -g según el tipo, pero aquí reforzamos por si acaso.
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
  target_compile_options(main PRIVATE
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:RelWithDebInfo>:-O2;-g>
  )
endif()

# Mensajes útiles
message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
message(STATUS "Usando TBB target: ${_TBB_TARGET}")
