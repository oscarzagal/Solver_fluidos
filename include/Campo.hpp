//
// Created by oscar on 22/06/25.
//

#ifndef CAMPOS_HPP
#define CAMPOS_HPP

#include <vector>

#include "condiciones_de_frontera.hpp"
#include "condiciones_de_frontera_MDot.hpp"

namespace Campo {

struct velFace {
    // Velocidades en las caras en "x" e "y"
    std::vector<double> uFace_x, vFace_y;

    // Constructor
    velFace(int, int);
};

struct MDotStar {
    // Flujos de masa en "x" e "y"
    std::vector<double> mDotStar_x, mDotStar_x_old;
    std::vector<double> mDotStar_y, mDotStar_y_old;

    // Listas de condiciones de frontera
    std::vector<CF_MDot<Dirichlet_MDot>>    lista_Dirichlet_x;
    std::vector<CF_MDot<Zero_Neumann_MDot>> lista_Zero_Neumann_x;

    std::vector<CF_MDot<Dirichlet_MDot>>    lista_Dirichlet_y;
    std::vector<CF_MDot<Zero_Neumann_MDot>> lista_Zero_Neumann_y;

    // Constructor
    MDotStar(int, int);
};

struct A_coef {
    // Matriz A que almacena los coeficientes agrupados
    std::vector<double> ap, ae, aw, an, as, b;
};

struct Momentum {
    std::vector<double> u_star, u_old; // Velocidad en x
    std::vector<double> v_star, v_old; // Velocidad en y
    MDotStar mdotstar;
    A_coef A_u, A_v;

    // Listas de condiciones de frontera (PROVISIONAL: cambiaran los tipos luego de
    // modificar la interfaz de condiciones de frontera)
    std::vector<Condicion_frontera::Dirichlet>             lista_parches_dirichlet_u;
    std::vector<std::shared_ptr<Condicion_frontera::Base>> lista_parches_dinamicos_u;

    std::vector<Condicion_frontera::Dirichlet>             lista_parches_dirichlet_v;
    std::vector<std::shared_ptr<Condicion_frontera::Base>> lista_parches_dinamicos_v;

    // Constructor
    Momentum(int, int);

};

struct Presion {
    std::vector<double> P_star, P_old; // Presión
    A_coef A_p;

    // Listas de condiciones de frontera (PROVISIONAL: cambiaran los tipos luego de
    // modificar la interfaz de condiciones de frontera)
    std::vector<Condicion_frontera::Dirichlet>             lista_parches_dirichlet;
    std::vector<std::shared_ptr<Condicion_frontera::Base>> lista_parches_dinamicos;

    // Constructor
    Presion(int, int);
};

} // Fin namespace Campo_Escalar



#endif //CAMPOS_HPP
