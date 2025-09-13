//
// Created by oscar on 22/06/25.
//

#ifndef ECUACIONES_GOBERNANTES_HPP
#define ECUACIONES_GOBERNANTES_HPP

#include <vector>

namespace Ecuaciones_gobernantes {


struct MDotStar {
    // Flujos de masa en "x" e "y"
    std::vector<double> mDotStar_x, mDotStar_x_old;
    std::vector<double> mDotStar_y, mDotStar_y_old;
};

struct A_coef {
    // Matriz A que almacena los coeficientes agrupados
    std::vector<double> ac, ae, aw, an, as, b;
};

struct Momentum {
    std::vector<double> u_star, u_old; // Velocidad en x
    std::vector<double> v_star, v_old; // Velocidad en y
    MDotStar mdotstar;
    A_coef A_u, A_v;

    // Constructor
    Momentum(const int, const int);

};

struct Presion {
    std::vector<double> P_star, P_old; // Presión
    A_coef A_p;

    // Constructor
    Presion(const int, const int);
};


} // Fin namespace Ecuaciones_gobernantes

#endif //ECUACIONES_GOBERNANTES_HPP
