//
// Created by oscar on 11/09/25.
//

#include "ecuacion_momentum.hpp"
#include "esquemas_de_discretizacion.hpp"
#include "malla_por_bloques.hpp"
#include <iostream>
#include <vector>

Ecuacion_Momentum::Ecuacion_Momentum
(
    Malla::Mallador    & malla_,
    Campo::Momentum    & velU_,
    Campo::Presion     & presion_,
    Gradiente          & grad_,
    fluxes_difusivos   & flux_dif_,
    fluxes_convectivos & flux_con_v_
) :
    malla(malla_),
    velU(velU_),
    presion(presion_),
    grad(grad_),
    flux_dif(flux_dif_),
    flux_conv(flux_con_v_),
    nx(malla_.obtener_el_numero_de_nodos(Malla::Nodos::nx)),
    ny(malla_.obtener_el_numero_de_nodos(Malla::Nodos::ny)),
    inter(Malla::Mallador::obtener_factores_de_interpolacion(malla_))
{}

void Ecuacion_Momentum::resolver() {

    std::cout << "nx: " << nx << "\n";
    std::cout << "ny: " << ny << "\n";

    std::vector<double> vol = malla.obtener_volumenes();

    std::cout << "ge en el ecuacion_momentum.cpp \n";
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            printf("ge[%d] = %f\n", i + nx * j, inter.ge[i + nx * j]);
        }
    }

    std::cout << "\n\n";

    std::cout << "volumenes en ecuacion_momentum.cpp \n";
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            printf("vol[%d] = %f\n", i + nx * j, vol[i + nx * j]);
        }
    }

    Discretizacion::Explicita::gradiente(nx, ny, inter, grad, presion.P_star, malla);

}

