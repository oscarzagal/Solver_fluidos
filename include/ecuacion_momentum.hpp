//
// Created by oscar on 11/09/25.
//

#ifndef ECUACION_MOMENTUM_HPP
#define ECUACION_MOMENTUM_HPP

#include "Campo.hpp"
#include "variables_discretizacion.hpp"

struct Ecuacion_Momentum{

    // Constructor
    Ecuacion_Momentum
    (
        Malla::Mallador    &,
        Campo::Momentum    &,
        Campo::Presion     &,
        Gradiente          &,
        fluxes_difusivos   &,
        fluxes_convectivos &
    );

    /* Funciones miembro */

    // Modifica el estado de "velU"
    void resolver();


    /* Miembros */

    // Parametros del constructor
    Malla::Mallador    & malla;
    Campo::Momentum    & velU;
    Campo::Presion     & presion;
    Gradiente          & grad;
    fluxes_difusivos   & flux_dif;
    fluxes_convectivos & flux_conv;

    // Miembros adicionales que se inicializan en la lista de inicializacion
    const int nx;
    const int ny;
    const Malla::Mallador::Interpolacion inter;

};


#endif //ECUACION_MOMENTUM_HPP
