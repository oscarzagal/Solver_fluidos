//
// Created by oscar on 12/06/25.
//

#ifndef ESQUEMAS_DE_DISCRETIZACION_HPP
#define ESQUEMAS_DE_DISCRETIZACION_HPP

#include "malla_por_bloques.hpp"
#include "ecuaciones_gobernantes.hpp"
#include "variables_discretizacion.hpp"

namespace Esquemas_discretizacion {

    void laplaciano_lineal
    (
        int nx,
        int ny,
        double gamma,
        fluxes_difusivos &fluxes,
        const Malla::Mallador &malla
    );

    void gradiente_explicito
    (
        int nx,
        int ny,
        const Malla::Mallador::Interpolacion &inter,
        Grad_explicito &grad_explicito,
        const std::vector<double> &P,
        const Malla::Mallador &malla
    );

    void divergencia_upwind
    (
        int nx,
        int ny,
        fluxes_convectivos &fluxes,
        Ecuaciones_gobernantes::MDotStar &mstar
    );

    void divergencia_explicita
    (
        int nx,
        int ny,
        const Malla::Mallador::Interpolacion &inter,
        fluxes_convectivos &fluxes,
        const std::vector<double> &P,
        const Malla::Mallador &malla
    );

    void construccion_matriz_A_momentum
    (
        int nx,
        int ny,
        const fluxes_difusivos &fluxes_dif,
        const fluxes_convectivos &fluxes_conv,
        Ecuaciones_gobernantes::A_coef &A_u,
        Ecuaciones_gobernantes::A_coef &A_v,
        Grad_explicito &grad_explicito
    );

    void construccion_coeficiente_b_momemtum
    (
        int nx,
        int ny,
        Ecuaciones_gobernantes::A_coef &A_u,
        Ecuaciones_gobernantes::A_coef &A_v,
        const Malla::Mallador &malla,
        const Grad_explicito &grad_explicito,
        const double rho
    );


} // Fin namespace Esquemas_discretizacion


#endif //ESQUEMAS_DE_DISCRETIZACION_HPP
