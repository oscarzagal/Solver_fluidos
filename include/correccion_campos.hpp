//
// Created by oscar on 04/10/25.
//

#ifndef SOLVER_FLUIDOS_CORRECCION_CAMPOS_HPP
#define SOLVER_FLUIDOS_CORRECCION_CAMPOS_HPP

#include "Campo.hpp"
#include "flujo_de_masa.hpp"
#include "malla_por_bloques.hpp"
#include <vector>


struct Correccion {

    // Constructor
    Correccion
    (
        const Malla::Mallador &,
        const Coeficiente_d   &,
        Campo::MDotStar       &,
        Campo::Momentum       &,
        Campo::Presion        &
    );

    /* Funciones miembro */

    // Modifica a "mdotstar", "vel" y "presion"
    void corregir();

    // Modifica a "celdas_interiores". Esto es util para usar a "for_each"
    void obtener_celdas_interiores();


    /* Miembros */

    // Parametros del constructor
    const Malla::Mallador & malla;
    const Coeficiente_d   & coef_d;
    Campo::MDotStar       & mdotstar;
    Campo::Momentum       & velU;
    Campo::Presion        & presion;


    // Miembros adicionales que se inicializan en la lista de inicializacion
    const int nx;
    const int ny;
    const Malla::Mallador::Interpolacion inter;

    // Instancia de gradiente para el calculo de los gradientes de presion
    // corregida
    Gradiente gradPrime_times_vol;

    // Volumen de los elementos computacionales
    std::vector<double> vol;

    // Celdas computacionales
    std::vector<int> celdas_interiores;

};


#endif //SOLVER_FLUIDOS_CORRECCION_CAMPOS_HPP
