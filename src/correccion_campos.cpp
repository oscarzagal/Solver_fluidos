//
// Created by oscar on 04/10/25.
//

#include "correccion_campos.hpp"
#include "utilidades.hpp"
#include "variables_discretizacion.hpp"
#include "esquemas_de_discretizacion.hpp"
#include <algorithm>
#include <execution>
#include <iostream>

// Constructor
Correccion::Correccion
(
    const Malla::Mallador & malla_,
    const Coeficiente_d   & coef_d_,
    Campo::MDotStar       & mdotstar_,
    Campo::Momentum       & velU_,
    Campo::Presion        & presion_
 ) :
    malla(malla_),
    coef_d(coef_d_),
    mdotstar(mdotstar_),
    velU(velU_),
    presion(presion_),
    nx(malla_.obtener_el_numero_de_nodos(Malla::Nodos::nx)),
    ny(malla_.obtener_el_numero_de_nodos(Malla::Nodos::ny)),
    inter(Malla::Mallador::obtener_factores_de_interpolacion(malla_)),
    gradPrime_times_vol(nx, ny),
    vol(malla_.obtener_volumenes()),
    celdas_interiores((nx - 2) * (ny - 2), 0)

{}


/*-----------------------------------------------------------------------------
                            Funciones miembro
-----------------------------------------------------------------------------*/

void Correccion::obtener_celdas_interiores() {

    int iter = 0;
    for (int j = 1 ; j < ny - 1 ; ++j) {
      for (int i = 1 ; i < nx - 1 ; ++i) {
          celdas_interiores[iter] = idx(i, j, nx);
          iter++;
      }
    }

}


void Correccion::corregir() {

    Discretizacion::Explicita::gradiente(nx, ny, inter, gradPrime_times_vol, presion.Pprime, malla);

    auto& u_star = velU.u_star;
    auto& v_star = velU.v_star;

    const auto& dC_u = coef_d.dC_u;
    const auto& dC_v = coef_d.dC_v;

    // TODO: poner tambien dentro el calculo de la correccion del flujo de masa
    // Calculo de las velocidades corregidas
    auto bucle = [&](int Centro) {

        const double gradPrime_x = gradPrime_times_vol.grad_x_vol[Centro] / vol[Centro];
        const double gradPrime_y = gradPrime_times_vol.grad_y_vol[Centro] / vol[Centro];

        u_star[Centro] = u_star[Centro] - dC_u[Centro] * gradPrime_x;
        v_star[Centro] = v_star[Centro] - dC_v[Centro] * gradPrime_y;

    };

    std::cout << "celdas_interiores.begin() = " << *celdas_interiores.begin() << "\n";
    std::cout << "celdas_interiores.end() - 1 = " << *(celdas_interiores.end() - 1) << "\n";

    // Paralelizacion con "for_each"
    std::for_each
    (
        std::execution::par_unseq,
        celdas_interiores.begin(),
        celdas_interiores.end(),
        bucle
    );







}


/*-----------------------------------------------------------------------------
                            Fin Funciones miembro
-----------------------------------------------------------------------------*/
