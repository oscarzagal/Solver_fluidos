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

# Para ejecutar todos tests:
# ctest --test-dir cmake-build-debug/ --output-on-failure -j

cmake_minimum_required(VERSION 3.10)
project(solver_fluidos VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "No build type selected, defaulting to Debug")
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

# === Fuentes ===
# NO metas main.cpp en la librería.
file(GLOB_RECURSE SRC_FILES CONFIGURE_DEPENDS src/*.cpp)
list(REMOVE_ITEM SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)

add_library(core ${SRC_FILES})
target_include_directories(core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(main main.cpp)
target_link_libraries(main PRIVATE core)

# ===== oneTBB =====
find_package(TBB QUIET CONFIG)
if(NOT TBB_FOUND)
  find_package(TBB REQUIRED)
endif()

set(_TBB_TARGET "")
if(TARGET TBB::tbb)
  set(_TBB_TARGET TBB::tbb)
elseif(TARGET tbb)
  set(_TBB_TARGET tbb)
endif()

if(_TBB_TARGET STREQUAL "")
  message(FATAL_ERROR "No se encontró un target de TBB exportado. Revisa que 'libtbb-dev' esté instalado.")
endif()

# Enlaza TBB en core (PUBLIC para que main lo herede)
target_link_libraries(core PUBLIC ${_TBB_TARGET})

# ===== Definiciones y flags por configuración =====
foreach(tgt IN ITEMS core main)
  target_compile_definitions(${tgt} PRIVATE
    $<$<CONFIG:Debug>:DEBUG>
    $<$<CONFIG:Release>:NDEBUG;PARALLEL>
    $<$<CONFIG:RelWithDebInfo>:NDEBUG;PARALLEL>
  )
  target_compile_options(${tgt} PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall;-Wextra;-Wpedantic>
    $<$<AND:$<CONFIG:Debug>,$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>>:-fsanitize=address,undefined>
  )
  target_link_options(${tgt} PRIVATE
    $<$<AND:$<CONFIG:Debug>,$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>>:-fsanitize=address,undefined>
  )
endforeach()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
  target_compile_options(core PRIVATE
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:RelWithDebInfo>:-O2;-g>
  )
  target_compile_options(main PRIVATE
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:RelWithDebInfo>:-O2;-g>
  )
endif()

message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
message(STATUS "Usando TBB target: ${_TBB_TARGET}")

# === Tests ===
include(FetchContent)
FetchContent_Declare(
  catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.6.0
)
FetchContent_MakeAvailable(catch2)

enable_testing()
add_subdirectory(tests)
