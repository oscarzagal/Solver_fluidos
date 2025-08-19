//
// Created by oscar on 22/06/25.
//

#include "config_control.hpp"
#include "esquemas_de_discretizacion.hpp"
#include "malla_por_bloques.hpp"
#include "ecuaciones_gobernantes.hpp"
#include <vector>


namespace Ecuaciones_gobernantes {

/*-----------------------------------------------------------------------------
                            Ecuacion de Energia
-----------------------------------------------------------------------------*/

Energia::Energia(Malla::Mallador & malla_) : malla(malla_) {}

void Energia::ensamblar() {

    const int nx = static_cast<int>(malla.obtener_coordenadas_tmp_x().size());
    const int ny = static_cast<int>(malla.obtener_coordenadas_tmp_y().size());

    A_coef A;

    // NOTE: uso no correcto del esquema de discretizacion
    // Esquemas_discretizacion::laplaciano_lineal(nx,ny,k,A,malla);

    asignar_matriz(A);

}

void Energia::asignar_matriz(const A_coef& A_paso) {
    A = A_paso;
}


A_coef Energia::obtener_coeficientes() {
    return A;
}


/*-----------------------------------------------------------------------------
                            Ecuacion de Momentum
-----------------------------------------------------------------------------*/

Momentum::Momentum(const double nu_, Malla::Mallador & malla_, std::vector<double> & Pstar_) :

    // Obtencion de los nodos una unica vez en el tiempo de vida del objeto
    nx(malla_.obtener_el_numero_de_nodos(Malla::Nodos::nx)),
    ny(malla_.obtener_el_numero_de_nodos(Malla::Nodos::ny)),
    // Inicializacion del flujo de masa
    mstar{std::vector<double>(nx*ny,0.0), std::vector<double>(nx*ny,0.0)},
    nu(nu_),
    malla(malla_),
    inter(Malla::Mallador::obtener_factores_de_interpolacion(malla_)),
    Pstar(Pstar_)
{}

void Momentum::unir_ecuacion() {

    // Modifica el estado de "gradP_explicito"
    Esquemas_discretizacion::gradiente_explicito(nx,ny,inter,gradP_explicito,Pstar,malla);

    // Modifica el estado de "fluxes_conv"
    Esquemas_discretizacion::divergencia_upwind(nx,ny,flux_conv,mstar);

    // Modifica el estado de "flux_dif"
    Esquemas_discretizacion::laplaciano_lineal(nx,ny,nu,flux_dif,malla);

    // TODO: implementar la funcion de abajo
    // Modifica el estado de "A_u" y "A_v"
    Esquemas_discretizacion::construccion_matriz_A_momentum(nx,ny,flux_dif,flux_conv,A_u,A_v,gradP_explicito);

    // TODO: hacer una funcion para el calculo del flujo de masa

}

// TODO: hacer enums para estas madres
A_coef Momentum::obtener_coeficientes_para_u() const {
    return A_u;
}

A_coef Momentum::obtener_coeficientes_para_v() const {
    return A_v;
}

void Momentum::asignar_matriz_para_u(const A_coef& A_paso) {
    A_u = A_paso;
}

void Momentum::asignar_matriz_para_v(const A_coef& A_paso) {
    A_v = A_paso;
}


} // Fin namespace Ecuaciones_gobernantes
