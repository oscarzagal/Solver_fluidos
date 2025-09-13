//
// Created by oscar on 11/09/25.
//

#include "ecuacion_momentum.hpp"
#include "Campo.hpp"
#include "config_control.hpp"
#include "esquemas_de_discretizacion.hpp"
#include "malla_por_bloques.hpp"
#include "solvers_lineales.hpp"
#include <iostream>
#include <variant>

/*-----------------------------------------------------------------------------
                                 Constructor
-----------------------------------------------------------------------------*/

Ecuacion_Momentum::Ecuacion_Momentum
(
    Malla::Mallador    & malla_,
    Campo::Momentum    & velU_,
    Campo::Presion     & presion_,
    Campo::MDotStar    & mdotstar_,
    Gradiente          & grad_,
    fluxes_difusivos   & flux_dif_,
    fluxes_convectivos & flux_con_v_
) :
    malla(malla_),
    velU(velU_),
    presion(presion_),
    mdotstar(mdotstar_),
    grad(grad_),
    flux_dif(flux_dif_),
    flux_conv(flux_con_v_),
    nx(malla_.obtener_el_numero_de_nodos(Malla::Nodos::nx)),
    ny(malla_.obtener_el_numero_de_nodos(Malla::Nodos::ny)),
    vol(malla_.obtener_volumenes()),
    inter(Malla::Mallador::obtener_factores_de_interpolacion(malla_)),
    solver_u(Solver_lineal::solverElegido(nx, ny, lambda_Vel, velU.u_star, velU.u_old, solver_elegido_u)),
    solver_v(Solver_lineal::solverElegido(nx, ny, lambda_Vel, velU.v_star, velU.v_old, solver_elegido_v))
{}

/*-----------------------------------------------------------------------------
                             Fin Constructor
-----------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------
                             Funciones miembro
-----------------------------------------------------------------------------*/

void Ecuacion_Momentum::calcular_conductancia_difusiva(const double nu) {

    Discretizacion::Implicita::laplaciano_lineal(nx, ny, nu, flux_dif, malla);

}

void Ecuacion_Momentum::resolver() {

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


    // NOTE: Si la lambda está dentro de un método y solo accede a miembros: usa [this]
    // Solucion del sistema de ecuaciones
    std::visit([this](auto& s){
            s.resolver(velU.A_u);
            }, solver_u);

    std::visit([this](auto& s){
            s.resolver(velU.A_v);
            }, solver_v);


    for (int j = 0 ; j < ny ; ++j) {
      for (int i = 0 ; i < nx ; ++i) {
          printf("vel_u[%d] = %f\n", i + nx * j, velU.u_star[i + nx * j]);
      }
    }



    // TODO: Llamar al solver lineal

    // std::cout << "nx: " << nx << "\n";
    // std::cout << "ny: " << ny << "\n";

    // std::vector<double> vol = malla.obtener_volumenes();

    // std::cout << "ge en el ecuacion_momentum.cpp \n";
    // for (int j = 0; j < ny; ++j) {
    //     for (int i = 0; i < nx; ++i) {
    //         printf("ge[%d] = %f\n", i + nx * j, inter.ge[i + nx * j]);
    //     }
    // }

    // std::cout << "\n\n";

    // std::cout << "volumenes en ecuacion_momentum.cpp \n";
    // for (int j = 0; j < ny; ++j) {
    //     for (int i = 0; i < nx; ++i) {
    //         printf("vol[%d] = %f\n", i + nx * j, vol[i + nx * j]);
    //     }
    // }



}


/*-----------------------------------------------------------------------------
                         Fin Funciones miembro
-----------------------------------------------------------------------------*/
