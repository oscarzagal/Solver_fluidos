//
// Created by oscar on 28/09/25.
//

#ifndef ECUACION_PRESION_HPP
#define ECUACION_PRESION_HPP

#include "Campo.hpp"
#include "malla_por_bloques.hpp"
#include "flujo_de_masa.hpp"
#include "solvers_lineales.hpp"
#include "variables_discretizacion.hpp"
#include <vector>

struct Ecuacion_Presion {

    // Constructor
    Ecuacion_Presion
    (
        const Malla::Mallador &,
        Campo::Presion        &,
        const Campo::MDotStar &,
        const Coeficiente_d   &
    );


    /* Funciones miembro */

    // Resolucion del campo de presion de correccion. Modifica el estado de
    // "presion.Pprime"
    void resolver();


    /* Miembros */

    // Parametros del constructor
    const Malla::Mallador & malla;
    Campo::Presion        & presion;
    const Campo::MDotStar & mdotstar;
    const Coeficiente_d   & coef_d;

    // Miembros adicionales que se inicializan en la lista de inicializacion
    const int nx;
    const int ny;
    const std::vector<double> vol; // Volumenes de las celdas

    // Eleccion del solver lineal
    Solver_lineal::solverVariant solver_Pprime;

    // Fluxes difusivos para la ecuacion de presion
    fluxes_difusivos flux_dif_P;

    // Factores de interpolacion
    const Malla::Mallador::Interpolacion inter;

};


#endif //ECUACION_PRESION_HPP
