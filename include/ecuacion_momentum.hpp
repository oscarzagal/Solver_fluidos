//
// Created by oscar on 11/09/25.
//

#ifndef ECUACION_MOMENTUM_HPP
#define ECUACION_MOMENTUM_HPP

#include "Campo.hpp"
#include "variables_discretizacion.hpp"
#include "solvers_lineales.hpp"
#include "flujo_de_masa.hpp"
#include <vector>

struct Ecuacion_Momentum {

    // Constructor
    Ecuacion_Momentum
    (
        Malla::Mallador    &,
        Campo::Momentum    &,
        Campo::Presion     &,
        Campo::MDotStar    &,
        Gradiente          &,
        fluxes_difusivos   &,
        fluxes_convectivos &
    );

    /* Funciones miembro */

    // Modifica el estado de "velU"
    void resolver();

    // Modifica el estado de "flux_dif". Para flujo incompresible solo es
    // necesario calcularla una vez.
    void calcular_conductancia_difusiva();


    // Helper para el solver lineal
    template<typename Variant, class Mat>
        void resolver_con(Variant& var, Mat& A) {
            std::visit([&](auto& campo){ campo.resolver(A); }, var);
        }


    /* Miembros */

    // Parametros del constructor
    Malla::Mallador    & malla;
    Campo::Momentum    & velU;
    Campo::Presion     & presion;
    Campo::MDotStar    & mdotstar;
    Gradiente          & grad;
    fluxes_difusivos   & flux_dif;
    fluxes_convectivos & flux_conv;

    // Miembros adicionales que se inicializan en la lista de inicializacion
    const int nx;
    const int ny;
    const std::vector<double> vol; // Volumenes de las celdas
    const Malla::Mallador::Interpolacion inter;

    // Eleccion del solver lineal
    Solver_lineal::solverVariant solver_u;
    Solver_lineal::solverVariant solver_v;

    // Coeficiente "d"
    Coeficiente_d coef_d;

};


#endif //ECUACION_MOMENTUM_HPP
