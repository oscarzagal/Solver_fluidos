//
// Created by oscar on 12/06/25.
//

#ifndef ESQUEMAS_DE_DISCRETIZACION_HPP
#define ESQUEMAS_DE_DISCRETIZACION_HPP

#include "malla_por_bloques.hpp"
#include "variables_discretizacion.hpp"
#include "Campo.hpp"
#include "flujo_de_masa.hpp"
#include <vector>

namespace Discretizacion {

    namespace Implicita {

        // Modifica el estado de "fluxes"
        void laplaciano_lineal
        (
             int nx,
             int ny,
             double gamma,
             fluxes_difusivos      & fluxes,
             const Malla::Mallador & malla
        );

        // Modifica el estado de "fluxes"
        void laplaciano_lineal_presion
        (
            int nx,
            int ny,
            const Coeficiente_d                  & coef_d,
            fluxes_difusivos                     & fluxes,
            const Malla::Mallador::Interpolacion & inter,
            const Malla::Mallador                & malla
        );

        // Modifica el estado de "fluxes"
        void divergencia_upwind
        (
             int nx,
             int ny,
             fluxes_convectivos & fluxes,
             Campo::MDotStar    & mdotstar
        );
    }

    namespace Explicita {

        // Modifica el estado de "grad_explicito". Esta funcion calcula el gradiente
        // por el volumen del elemento computacional.
        void gradiente
        (
             int nx,
             int ny,
             const Malla::Mallador::Interpolacion &inter,
             Gradiente &grad_explicito,
             const std::vector<double> &P,
             const Malla::Mallador &malla
        );

    }

    // Modifica el estado de "A", espeficamente al miembro "b"
    void construccion_coeficiente_b_momemtum
    (
        int nx,
        int ny,
        Campo::A_coef & A,
        const std::vector<double>      & vol,
        const std::vector<double>      & grad_vol,
        const std::vector<double>      & vel
    );

    // Modifica a "A_u" y "A_v"
    void construccion_matriz_A_momentum
    (
        int nx,
        int ny,
        const fluxes_difusivos    & fluxes_dif,
        const fluxes_convectivos  & fluxes_conv,
        const std::vector<double> & vol,
        const std::vector<double> & vel_u,
        const std::vector<double> & vel_v,
        Campo::A_coef             & A_u,
        Campo::A_coef             & A_v,
        Gradiente                 & grad
    );

    // Modifica a "A_p"
    void construccion_matriz_A_presion
    (
        int nx,
        int ny,
        const fluxes_difusivos & fluxes_difusivos,
        const Campo::MDotStar  & mdotstar,
        Campo::A_coef          & A_p
    );


} // Fin namespace Esquemas_discretizacion


#endif //ESQUEMAS_DE_DISCRETIZACION_HPP
