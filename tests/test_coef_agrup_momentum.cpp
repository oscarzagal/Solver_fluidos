#include "catch2/catch_approx.hpp"
#include "catch2/catch_test_macros.hpp"
#include <catch2/catch_all.hpp>
#include <vector>
#include "esquemas_de_discretizacion.hpp"
#include "malla_por_bloques.hpp"
#include "variables_discretizacion.hpp"

TEST_CASE("coeficientes agrupados momentum", "[coeficientes_momentum]") {

    // Variables globales:
    // inline constexpr double rho = 1.0;
    // inline constexpr double nu = 2.5e-3;
    // inline constexpr double delta_t = 0.5;
    // inline constexpr double lambda_Vel = 0.6;

    const int nx = 5, ny = 5;
    const double gamma = 2.5e-3;
    fluxes_difusivos flux_dif(nx, ny);
    fluxes_convectivos flux_conv(nx, ny);
    Gradiente grad(nx, ny);
    Campo::Presion presion(nx, ny, 1.0);
    Campo::MDotStar mdotstar(nx, ny, 0.0);
    Campo::Momentum velU(nx, ny, 0.0005, 0.0);

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

    const Malla::Mallador::Interpolacion inter = Malla::Mallador::obtener_factores_de_interpolacion(malla);

    std::vector vol = malla.obtener_volumenes();

    Discretizacion::Implicita::laplaciano_lineal(nx, ny, gamma, flux_dif, malla);

    Discretizacion::Explicita::gradiente(nx, ny, inter, grad, presion.P_star, malla);

    Discretizacion::Implicita::divergencia_upwind(nx, ny, flux_conv, mdotstar);

    Discretizacion::construccion_matriz_A_momentum
    (
        nx,
        ny,
        flux_dif,
        flux_conv,
        vol,
        velU.u_star,
        velU.v_star,
        velU.A_u,
        velU.A_v,
        grad
    );

    std::vector<double> ac_u = {1.0000000000000000,  1.0000000000000000,  1.0000000000000000,  1.0000000000000000,  1.0000000000000000,  1.0000000000000000, 0.39537037037037037, 0.39120370370370372, 0.39537037037037037,  1.0000000000000000,  1.0000000000000000, 0.39120370370370372, 0.38703703703703707, 0.39120370370370372,  1.0000000000000000,  1.0000000000000000, 0.39537037037037032, 0.39120370370370372, 0.39537037037037032,  1.0000000000000000,  1.0000000000000000,  1.0000000000000000,  1.0000000000000000,  1.0000000000000000,  1.0000000000000000};
    std::vector<double> b_u = {0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 1.9018518518518516E-004, 1.8935185185185184E-004, 1.9018518518518516E-004, 0.0000000000000000     , 0.0000000000000000     , 1.8935185185185184E-004, 1.8851851851851853E-004, 1.8935185185185184E-004, 0.0000000000000000     , 0.0000000000000000     , 1.9018518518518516E-004, 1.8935185185185184E-004, 1.9018518518518516E-004, 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000};
    std::vector<double> ac_v = {1.0000000000000000,   1.0000000000000000,   1.0000000000000000,   1.0000000000000000,   1.0000000000000000,   1.0000000000000000,  0.39537037037037037,  0.39120370370370372,  0.39537037037037037,   1.0000000000000000,   1.0000000000000000,  0.39120370370370372,  0.38703703703703707,  0.39120370370370372,   1.0000000000000000,   1.0000000000000000,  0.39537037037037032,  0.39120370370370372,  0.39537037037037032,   1.0000000000000000,   1.0000000000000000,   1.0000000000000000,   1.0000000000000000,   1.0000000000000000,   1.0000000000000000};
    std::vector<double> b_v = {0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000};

    REQUIRE(velU.A_u.ac.size() == static_cast<int>(nx * ny));
    REQUIRE(velU.A_u.b.size()  == static_cast<int>(nx * ny));
    REQUIRE(velU.A_v.ac.size() == static_cast<int>(nx * ny));
    REQUIRE(velU.A_v.b.size()  == static_cast<int>(nx * ny));

    for (int j = 0 ; j < ny ; ++j) {
        for (int i = 0 ; i < nx ; ++i) {

            const int Centro = i + nx * j;

            REQUIRE(velU.A_u.ac[Centro] == Catch::Approx(ac_u[Centro]));
            REQUIRE(velU.A_u.b[Centro]  == Catch::Approx(b_u[Centro]));
            REQUIRE(velU.A_v.ac[Centro] == Catch::Approx(ac_v[Centro]));
            REQUIRE(velU.A_v.b[Centro]  == Catch::Approx(b_v[Centro]));

        }
    }


}

