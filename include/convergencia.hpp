//
// Created by oscar on 08/10/25.
//

#ifndef SOLVER_FLUIDOS_CONVERGENCIA_HPP
#define SOLVER_FLUIDOS_CONVERGENCIA_HPP

#include <utility>
#include <vector>
#include "Campo.hpp"

enum class Residual {
    presion,
    u,
    v,
    mdotstar_x,
    mdotstar_y
};

// Devuelve el error mayor de todos los campos y el campo en cuestion
std::pair<double, std::string> error_mayor
(
    int nx,
    int ny,
    const Campo::Momentum & velU,
    const Campo::Presion  & presion,
    const Campo::MDotStar & mdotstar,
    std::vector<double>   & error_mayor_por_campo
);

// Modifica el valor de "error_mayor_por_campo"
void error_por_campo
(
    int nx,
    int ny,
    Residual residual,
    const std::vector<double> & phi_new,
    const std::vector<double> & phi_old,
    std::vector<double>       & error_mayor_por_campo
);

#endif //SOLVER_FLUIDOS_CONVERGENCIA_HPP
