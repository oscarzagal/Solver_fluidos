#include "catch2/catch_approx.hpp"
#include "catch2/catch_test_macros.hpp"
#include <catch2/catch_all.hpp>
#include <variant>
#include <vector>
#include <iostream>
#include "esquemas_de_discretizacion.hpp"
#include "malla_por_bloques.hpp"
#include "variables_discretizacion.hpp"
#include "solvers_lineales.hpp"
#include "config_control.hpp"
#include "flujo_de_masa.hpp"

TEST_CASE("flujo de masa", "[flujo_de_masa]") {

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

    // Solo se debe usar una vez para evitar duplicaciones y cosas raras
    malla.preparar_parches_fronteras();

    typedef std::vector<Malla::Mallador::Parche> almacenar;

    // Obtener parches para su uso futuro
    almacenar Parches_norte = malla.obtener_parches(Malla::Frontera::Norte);
    almacenar Parches_sur   = malla.obtener_parches(Malla::Frontera::Sur);
    almacenar Parches_este  = malla.obtener_parches(Malla::Frontera::Este);
    almacenar Parches_oeste = malla.obtener_parches(Malla::Frontera::Oeste);

    construir_condiciones_de_frontera
    (
        Parches_norte,
        Parches_sur,
        Parches_este,
        Parches_oeste,
        velU.u_star,
        nx,
        velU.lista_parches_dirichlet_u,
        velU.lista_parches_dinamicos_u,
        g_dirichlet_u,
        g_zero_neumann_u
    );

    construir_condiciones_de_frontera
    (
        Parches_norte,
        Parches_sur,
        Parches_este,
        Parches_oeste,
        velU.v_star,
        nx,
        velU.lista_parches_dirichlet_v,
        velU.lista_parches_dinamicos_v,
        g_dirichlet_v,
        g_zero_neumann_v
    );

    // Aplicacion de las condiciones de frontera para u
    for (int i = 0 ; i < static_cast<int>(velU.lista_parches_dirichlet_u.size()) ; ++i) {
        velU.lista_parches_dirichlet_u[i].aplicar();
    }

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


    // Resolucion de los sistemas de ecuaciones para la velocidad
    Solver_lineal::solverVariant solver_u = Solver_lineal::solverElegido(nx, ny, lambda_Vel_SL, velU.u_star, velU.u_old, solver_elegido_u);
    Solver_lineal::solverVariant solver_v = Solver_lineal::solverElegido(nx, ny, lambda_Vel_SL, velU.v_star, velU.v_old, solver_elegido_v);

    std::visit([&](auto& campo_u) {
        campo_u.resolver(velU.A_u);
    }, solver_u);

    std::visit([&](auto& campo_v) {
        campo_v.resolver(velU.A_v);
    }, solver_v);


    // for (int j = 1 ; j < ny - 1 ; ++j) {
    //     for (int i = 1 ; i < nx - 1 ; ++i) {

    //         const int Centro = i + nx * j;
    //         const int Este   = (i + 1) + nx * j;
    //         const int Oeste  = (i - 1) + nx * j;
    //         const int Norte  = i + nx * (j + 1);
    //         const int Sur    = i + nx * (j - 1);

    //         velU.u_star[i+nx*j]=(velU.u_star[i+1+nx*j]*(-velU.A_u.ae[i+nx*j])+velU.u_star[i-1+nx*j]*(-velU.A_u.aw[i+nx*j])
    //             +velU.u_star[i+nx*(j+1)]*(-velU.A_u.an[i+nx*j])+velU.u_star[i+nx*(j-1)]*(-velU.A_u.as[i+nx*j])+velU.A_u.b[i+nx*j])
    //             /velU.A_u.ac[i+nx*j];
    //     }
    // }

    // Coeficiente "d"
    Coeficiente_d coef_d(nx, ny, vol);

    // Velocidad en las caras
    Campo::velFace velface(nx, ny, 0.0);

    calcular_flujo_de_masa
    (
        nx,
        ny,
        velface,
        mdotstar,
        coef_d,
        malla,
        velU,
        presion,
        grad,
        vol,
        inter,
        velU.A_u,
        velU.A_v
    );

    // Para el campo "u"
    for (int i = 0; i < static_cast<int>(velU.lista_parches_dinamicos_u.size()); ++i) {
        velU.lista_parches_dinamicos_u[i]->aplicar();
    }

    // Valores esperados
    std::vector<double> ac_u = {1.0000000000000000,  1.0000000000000000,  1.0000000000000000,  1.0000000000000000,  1.0000000000000000,  1.0000000000000000, 0.39537037037037037, 0.39120370370370372, 0.39537037037037037,  1.0000000000000000,  1.0000000000000000, 0.39120370370370372, 0.38703703703703707, 0.39120370370370372,  1.0000000000000000,  1.0000000000000000, 0.39537037037037032, 0.39120370370370372, 0.39537037037037032,  1.0000000000000000,  1.0000000000000000,  1.0000000000000000,  1.0000000000000000,  1.0000000000000000,  1.0000000000000000};
    std::vector<double> b_u = {0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 1.9018518518518516E-004, 1.8935185185185184E-004, 1.9018518518518516E-004, 0.0000000000000000     , 0.0000000000000000     , 1.8935185185185184E-004, 1.8851851851851853E-004, 1.8935185185185184E-004, 0.0000000000000000     , 0.0000000000000000     , 1.9018518518518516E-004, 1.8935185185185184E-004, 1.9018518518518516E-004, 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000};
    std::vector<double> u_star = {0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 4.8735362997658077E-004, 4.9352865041641834E-004, 4.8731271043588834E-004, 0.0000000000000000     , 0.0000000000000000     , 4.9352865041641834E-004, 4.9991639883313068E-004, 4.9352785466278424E-004, 0.0000000000000000     , 0.0000000000000000     , 1.3133682733855088E-002, 1.3355409982623091E-002, 1.3214969864544896E-002, 0.0000000000000000     , 1.0000000000000000     , 1.0000000000000000     , 1.0000000000000000     , 1.0000000000000000     , 1.0000000000000000};
    std::vector<double> me_star = {0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 1.6348038006549985E-004, 1.6347356014205109E-004, 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 1.6557417487492483E-004, 1.6557404224931913E-004, 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 4.4148487860796969E-003, 4.4283966411946644E-003, 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000     , 0.0000000000000000};


    // Camprobacion de tamaños
    REQUIRE(mdotstar.mDotStar_x.size() == static_cast<int>(nx * ny));
    REQUIRE(velU.u_star.size() == static_cast<int>(nx * ny));

    // Coeficientes agrupados
    for (int j = 0 ; j < ny ; ++j) {
        for (int i = 0 ; i < nx ; ++i) {

            const int Centro = i + nx * j;

            REQUIRE(velU.A_u.ac[Centro] == Catch::Approx(ac_u[Centro]));
            REQUIRE(velU.A_u.b[Centro]  == Catch::Approx(b_u[Centro]));

        }
    }

    // REQUIRE(velU.u_star[6] == Catch::Approx(u_star[6]));

    // Velocidades en u
    for (int j = 1 ; j < ny - 1; ++j) {
        for (int i = 1 ; i < nx - 1; ++i) {

            const int Centro = i + nx * j;

            std::cout << "i = " << i << ", " << "j = " << j << "\n";
            REQUIRE(velU.u_star[Centro] == Catch::Approx(u_star[Centro]));

        }
    }


    // Flujo de masa en la cara "e"
    for (int j = 1 ; j < ny - 1; ++j) {
        for (int i = 1 ; i < nx - 2; ++i) {

            const int Centro = i + nx * j;
            const int Este   = (i + 1) + nx * j;

            std::cout << "i = " << i << ", " << "j = " << j << "\n";
            REQUIRE(mdotstar.mDotStar_x[Este] == Catch::Approx(me_star[Centro]));

        }
    }

}



