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

/* Miembros */

std::vector<double> dC_u, dE_u, dW_u, dC_v, dN_v, dS_v;

// Volumen de los elementos computacionales
const std::vector<double>& vol;

const int nx, ny;


/* Funciones miembro */

// Calculo de los coeificientes "d_C = V_{C}/a_{C}", por ejemplo. Se usa dentro
// de "calcular_flujo_de_masa". Modifica a "coef_d".
void calcular(const Campo::A_coef& A_u, const Campo::A_coef& A_v);

Coeficiente_d(int, int, const std::vector<double>&);

};

// TODO: implementar funcion

// Modifica a "coef_d", "velFace", y "mdotstar". Calcula el flujo de masa a traves de la
// interpolacion de Rhie-Chow
void calcular_flujo_de_masa
(
     int nx,
     int ny,
     Campo::velFace                       & velFace,
     Campo::MDotStar                      & mdotstar,
     Coeficiente_d                        & coef_d,
     const Campo::Momentum                & vel,
     const Campo::Presion                 & Pstar,
     const Gradiente                      & gradP,
     const std::vector<double>            & vol,
     const Malla::Mallador::Interpolacion & inter,
     const Campo::A_coef                  & A_u,
     const Campo::A_coef                  & A_v
);


#endif //FLUJO_DE_MASA_HPP
