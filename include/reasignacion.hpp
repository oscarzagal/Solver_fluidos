//
// Created by oscar on 13/10/25.
//

#ifndef SOLVER_FLUIDOS_REASIGNACION_HPP
#define SOLVER_FLUIDOS_REASIGNACION_HPP

#include "Campo.hpp"

void reasignar
(
    Campo::Presion  & presion,
    Campo::Momentum & velU,
    Campo::MDotStar & mdotstar,
    Campo::velFace  & velface
);


#endif //SOLVER_FLUIDOS_REASIGNACION_HPP
