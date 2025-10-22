#include <catch2/catch_all.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <utility>
#include <vector>
#include "esquemas_de_discretizacion.hpp"
#include "malla_por_bloques.hpp"
#include "variables_discretizacion.hpp"

TEST_CASE("laplaciano lineal", "[laplaciano]") {

    const int nx = 5, ny = 5;
    const double gamma = 2.5e-3;
    fluxes_difusivos fluxes(nx, ny);

    std::vector<int> nodos_en_x = {nx};
    std::vector<int> nodos_en_y = {ny};

    std::vector<std::pair<double, double>> coordenadas_en_x = {
        {0.0, 1.0}
    };

    std::vector<std::pair<double, double>> coordenadas_en_y = {
        {0.0, 1.0}
    };

    std::vector<std::string> nombres_frontera_norte = {
        {"norte"},
    };

    std::vector<std::string> nombres_frontera_sur = {
        {"sur"},
    };

    std::vector<std::string> nombres_frontera_este = {
        {"este"},
    };

    std::vector<std::string> nombres_frontera_oeste = {
        {"oeste"},
    };

    Malla::Mallador malla
    (
        nodos_en_x,
        nodos_en_y,
        coordenadas_en_x,
        coordenadas_en_y,
        nombres_frontera_norte,
        nombres_frontera_sur,
        nombres_frontera_este,
        nombres_frontera_oeste
    );

    // Primero tienes que preparar para luego obtener las coordenadas persistentes
    malla.preparar_coordenadas_persistentes();


    // FUNCION A PROBAR
    Discretizacion::Implicita::laplaciano_lineal(nx, ny, gamma, fluxes, malla);

    // Donde se guarda la salida
    std::vector<double> fluxFDif_e = { 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 2.4999999999999996E-003, 2.5000000000000005E-003, 4.9999999999999975E-003, 0.0000000000000000     , 0.0000000000000000     , 2.4999999999999996E-003, 2.5000000000000005E-003, 4.9999999999999975E-003, 0.0000000000000000     , 0.0000000000000000     , 2.4999999999999996E-003, 2.5000000000000005E-003, 4.9999999999999975E-003, 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000};

    // lectura_archivo("~/Escritorio/Maestria/solver_fluidos/cmake-build-debug/De.dat", x, y, fluxFDif_e);

    REQUIRE(fluxes.fluxFDif_e.size() == static_cast<size_t>(nx * ny));
    REQUIRE(fluxFDif_e.size()        == static_cast<size_t>(nx * ny));

    for (int j = 1 ; j < ny - 1; ++j) {
      for (int i = 1 ; i < nx - 1; ++i) {

            const int Centro = i + nx * j;

            REQUIRE(fluxes.fluxFDif_e[Centro] == Catch::Approx(- fluxFDif_e[Centro]));

      }
    }


}

