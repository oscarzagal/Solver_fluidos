# Un binario de tests (puedes crear más: test_mdot.cpp, etc.)
add_executable(tests test_laplaciano_lineal.cpp)

target_link_libraries(tests
  PRIVATE
    core
    Catch2::Catch2WithMain
)

# Sanitizers solo en tests (útil)
target_compile_options(tests PRIVATE
  $<$<AND:$<CONFIG:Debug>,$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>>:-fsanitize=address,undefined -O1 -g>
)
target_link_options(tests PRIVATE
  $<$<AND:$<CONFIG:Debug>,$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>>>:-fsanitize=address,undefined>
)

# Descubre tests de Catch2 para ctest
include(Catch)
catch_discover_tests(tests)
