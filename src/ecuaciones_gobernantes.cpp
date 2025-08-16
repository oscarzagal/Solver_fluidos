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

    nu(nu_),
    // Obtencion de los nodos una unica vez en el tiempo de vida del objeto
    nx(malla_.obtener_el_numero_de_nodos(Malla::Nodos::nx)),
    ny(malla_.obtener_el_numero_de_nodos(Malla::Nodos::ny)),
    malla(malla_),
    Pstar(Pstar_)
{
    inter = Malla::Mallador::obtener_factores_de_interpolacion(malla_);
    // Inicializacion del flujo de masa
    mstar.mstar_x.resize(nx*ny,0.0);
    mstar.mstar_y.resize(nx*ny,0.0);
}

void Momentum::ensamblar() {

    // Modifica el estado de "gradP_explicito"
    Esquemas_discretizacion::gradiente_explicito(nx,ny,inter,gradP_explicito,Pstar,malla);

    // Modifica el estado de "fluxes_conv"
    Esquemas_discretizacion::divergencia_upwind(nx,ny,flux_conv,mstar);

    // Modifica el estado de "flux_dif"
    Esquemas_discretizacion::laplaciano_lineal(nx,ny,nu,flux_dif,malla);

    // Modifica el estado de "A"
    Esquemas_discretizacion::construccion_matriz_A(nx,ny,flux_dif,flux_conv,A);

    // TODO: obtener los coeficientes agrupados a traves de una funcion de union

}

void Momentum::asignar_matriz(const A_coef& A_paso) {
    A = A_paso;
}

const A_coef Momentum::obtener_coeficientes() const {
    return A;
}



} // Fin namespace Ecuaciones_gobernantes
