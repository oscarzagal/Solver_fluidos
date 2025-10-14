//
// Created by oscar on 11/09/25.
//

#include "ecuacion_momentum.hpp"
#include "Campo.hpp"
#include "config_control.hpp"
#include "esquemas_de_discretizacion.hpp"
#include "flujo_de_masa.hpp"
#include "malla_por_bloques.hpp"
#include "solvers_lineales.hpp"

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
    solver_u(Solver_lineal::solverElegido(nx, ny, lambda_Vel_SL, velU.u_star, velU.u_old, solver_elegido_u)),
    solver_v(Solver_lineal::solverElegido(nx, ny, lambda_Vel_SL, velU.v_star, velU.v_old, solver_elegido_v)),
    coef_d(nx, ny, vol),
    velface(nx, ny, 0.0)
{}

/*-----------------------------------------------------------------------------
                             Fin Constructor
-----------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------
                             Funciones miembro
-----------------------------------------------------------------------------*/

void Ecuacion_Momentum::calcular_conductancia_difusiva() {

    Discretizacion::Implicita::laplaciano_lineal(nx, ny, nu, flux_dif, malla);

}


void Ecuacion_Momentum::resolver() {

    /*-----------------------------------------------------------------------------
                            Resolucion ecuacion Momentum
    -----------------------------------------------------------------------------*/

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


    resolver_con(solver_u, velU.A_u);
    resolver_con(solver_v, velU.A_v);

    /*-----------------------------------------------------------------------------
                          Fin Resolucion ecuacion Momentum
    -----------------------------------------------------------------------------*/



    /*-----------------------------------------------------------------------------
                           Actulizacion flujo de masa
    -----------------------------------------------------------------------------*/

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


    /*-----------------------------------------------------------------------------
                         Fin Actulizacion flujo de masa
    -----------------------------------------------------------------------------*/



    /*-----------------------------------------------------------------------------
                        Actualizar condiciones de frontera
    -----------------------------------------------------------------------------*/

    // Para el campo "u"
    for (int i = 0; i < static_cast<int>(velU.lista_parches_dinamicos_u.size()); ++i) {
        velU.lista_parches_dinamicos_u[i]->aplicar();
    }

    // Para el campo "v"
    for (int i = 0; i < static_cast<int>(velU.lista_parches_dinamicos_v.size()); ++i) {
        velU.lista_parches_dinamicos_v[i]->aplicar();
    }

    /*-----------------------------------------------------------------------------
                        Actualizar condiciones de frontera
    -----------------------------------------------------------------------------*/

}



/*-----------------------------------------------------------------------------
                         Fin Funciones miembro
-----------------------------------------------------------------------------*/

