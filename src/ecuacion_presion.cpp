//
// Created by oscar on 28/09/25.
//

#include "ecuacion_presion.hpp"
#include "config_control.hpp"

/*-----------------------------------------------------------------------------
                               Ecuacion presion
-----------------------------------------------------------------------------*/

Ecuacion_Presion::Ecuacion_Presion
(
    Malla::Mallador & malla_,
    Campo::Presion  & presion_,
    Campo::MDotStar & mdotstar_,
    Coeficiente_d   & coef_d_
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
    gPprime_y(nx * ny, 0.0)
{}


/*-----------------------------------------------------------------------------
                            Fin Ecuacion presion
-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
                               Funciones miembro
-----------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------
                             Fin Funciones miembro
-----------------------------------------------------------------------------*/
