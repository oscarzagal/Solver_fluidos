//
// Created by oscar on 28/09/25.
//

#include "ecuacion_presion.hpp"
#include "config_control.hpp"
#include "variables_discretizacion.hpp"
#include "esquemas_de_discretizacion.hpp"

/*-----------------------------------------------------------------------------
                                 Constructor
-----------------------------------------------------------------------------*/

Ecuacion_Presion::Ecuacion_Presion
(
    const Malla::Mallador & malla_,
    Campo::Presion        & presion_,
    const Campo::MDotStar & mdotstar_,
    const Coeficiente_d   & coef_d_
) :
    malla(malla_),
    presion(presion_),
    mdotstar(mdotstar_),
    coef_d(coef_d_),
    nx(malla_.obtener_el_numero_de_nodos(Malla::Nodos::nx)),
    ny(malla_.obtener_el_numero_de_nodos(Malla::Nodos::ny)),
    vol(malla_.obtener_volumenes()),
    solver_Pprime(Solver_lineal::solverElegido(nx, ny, lambda_P_SL, presion.Pprime, presion.P_old, solver_elegido_P)),
    gPprime_x(nx * ny, 0.0),
    gPprime_y(nx * ny, 0.0),
    flux_dif_P(nx, ny),
    inter(Malla::Mallador::obtener_factores_de_interpolacion(malla_))
{}


/*-----------------------------------------------------------------------------
                             Fin Constructor
-----------------------------------------------------------------------------*/



/*-----------------------------------------------------------------------------
                               Funciones miembro
-----------------------------------------------------------------------------*/

void Ecuacion_Presion::resolver() {

    Discretizacion::Implicita::laplaciano_lineal_presion(nx, ny, coef_d, flux_dif_P, inter, malla);

    Discretizacion::construccion_matriz_A_presion(nx, ny, flux_dif_P, mdotstar, presion.A_p);

    // Resolucion de la ecuacion de correccione de presion
    std::visit([this](auto& campo_Pprime){
            campo_Pprime.resolver(this->presion.A_p);
            }, solver_Pprime);

}


/*-----------------------------------------------------------------------------
                             Fin Funciones miembro
-----------------------------------------------------------------------------*/
