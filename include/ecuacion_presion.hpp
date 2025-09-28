//
// Created by oscar on 28/09/25.
//

#ifndef ECUACION_PRESION_HPP
#define ECUACION_PRESION_HPP

#include "Campo.hpp"
#include "malla_por_bloques.hpp"
#include "flujo_de_masa.hpp"
#include "solvers_lineales.hpp"
#include <vector>

struct Ecuacion_Presion {

    // Constructor
    Ecuacion_Presion
    (
        Malla::Mallador &,
        Campo::Presion  &,
        Campo::MDotStar &,
        Coeficiente_d   &
    );


    /* Funciones miembro */

    // Resolucion del campo de presion de correccion. Modifica el estado de
    // "presion.Pprime"
    void resolver();

    // Calculo del gradiente de presion de correccion. Modifica el estado de
    // "gPprime_x" y "gPprime_y"
    void gradiente_Pprime();


    /* Miembros */

    // Parametros del constructor
    Malla::Mallador & malla;
    Campo::Presion  & presion;
    Campo::MDotStar & mdotstar;
    Coeficiente_d   & coef_d;

    // Miembros adicionales que se inicializan en la lista de inicializacion
    const int nx;
    const int ny;
    const std::vector<double> vol; // Volumenes de las celdas

    // Eleccion del solver lineal
    Solver_lineal::solverVariant solver_Pprime;

    // Gradiente de presion de correccion
    std::vector<double> gPprime_x;
    std::vector<double> gPprime_y;

};




#endif //ECUACION_PRESION_HPP
