//
// Created by oscar on 12/06/25.
//

#ifndef ESQUEMAS_DE_DISCRETIZACION_HPP
#define ESQUEMAS_DE_DISCRETIZACION_HPP

#include "malla_por_bloques.hpp"
#include "ecuaciones_gobernantes.hpp"

namespace Esquemas_discretizacion {

    struct grad_explicito {

        // Gradiente explicito multiplicado por el volumen
        std::vector<double> grad_x_vol, grad_y_vol;

    };

    struct fluxes_convectivos {
        // Coeficientes convectivos en las caras para el coeficiente agrupado F
        std::vector<double> fluxFConv_e,fluxFConv_w,fluxFConv_n,fluxFConv_s;

        // Coeficientes convectivos en las caras para el coeficiente agrupado C
        std::vector<double> fluxCConv_e,fluxCConv_w,fluxCConv_n,fluxCConv_s;
    };

    struct fluxes_difusivos {
        // Coeficientes difusivos en las caras para el coeficiente agrupado F
        std::vector<double> fluxFDif_e,fluxFDif_w,fluxFDif_n,fluxFDif_s;

        // Coeficientes difusivos en las caras para el coeficiente agrupado C
        std::vector<double> fluxCDif_e,fluxCDif_w,fluxCDif_n,fluxCDif_s;
    };

    void laplaciano_lineal
    (
        const int &nx,
        const int &ny,
        const double &k,
        Ecuaciones_gobernantes::A_coef &A,
        Malla::Mallador& malla
    );

    void gradiente_explicito
    (
        int nx,
        int ny,
        const Malla::Mallador::Interpolacion &inter,
        grad_explicito &grad_explicito,
        const std::vector<double> &P,
        const Malla::Mallador &malla
    );

    void divergencia_upwind
    (
        int nx,
        int ny,
        fluxes_convectivos &fluxes,
        const Ecuaciones_gobernantes::Momentum::mstar &mstar
    );

} // Fin namespace Esquemas_discretizacion


#endif //ESQUEMAS_DE_DISCRETIZACION_HPP