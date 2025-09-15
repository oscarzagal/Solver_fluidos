//
// Created by oscar on 10/09/25.
//

#ifndef FLUJO_DE_MASA_HPP
#define FLUJO_DE_MASA_HPP

#include "Campo.hpp"
#include "malla_por_bloques.hpp"
#include "variables_discretizacion.hpp"
#include <vector>

struct Coeficiente_d {

std::vector<double> dC_u, dE_u, dW_u, dC_v, dN_v, dS_v;

// TODO: implementar Constructor
Coeficiente_d(int, int);

};

// TODO: implementar funciones

// Calculo de los coeificientes "d_C = V_{C}/a_{C}", por ejemplo. Se usa dentro
// de "calcular_flujo_de_masa". Modifica a "coef_d".
void coeficiente_d
(
    Coeficiente_d &coef_d,
    const std::vector<double> &vol,
    const Campo::A_coef& A
);

// Modifica a "velFace", y "mdotstar". Calcula el flujo de masa a traves de la
// interpolacion de Rhie-Chow
void calcular_flujo_de_masa
(
    int nx,
    int ny,
    Campo::velFace                       & velFace,
    Campo::MDotStar                      & mdotstar,
    const Coeficiente_d                  & coef_d,
    const Campo::Momentum                & vel,
    const Campo::Presion                 & Pstar,
    const Gradiente                      & gradP,
    const std::vector<double>            & vol,
    const Malla::Mallador::Interpolacion & inter
);


#endif //FLUJO_DE_MASA_HPP
