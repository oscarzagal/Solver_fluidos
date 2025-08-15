//
// Created by oscar on 13/08/25.
//

#ifndef VARIABLES_DISCRETIZACION_HPP
#define VARIABLES_DISCRETIZACION_HPP

#include <vector>

struct Grad_explicito {
    // Gradiente explicito multiplicado por el volumen
    std::vector<double> grad_x_vol, grad_y_vol;
};

struct fluxes_convectivos {
    // Coeficientes convectivos en las caras para el coeficiente agrupado F
    std::vector<double> fluxFConv_e, fluxFConv_w, fluxFConv_n, fluxFConv_s;

    // Coeficientes convectivos en las caras para el coeficiente agrupado C
    std::vector<double> fluxCConv_e, fluxCConv_w, fluxCConv_n, fluxCConv_s;
};

struct fluxes_difusivos {
    // Coeficientes difusivos en las caras para el coeficiente agrupado F
    std::vector<double> fluxFDif_e, fluxFDif_w, fluxFDif_n, fluxFDif_s;

    // Coeficientes difusivos en las caras para el coeficiente agrupado C
    std::vector<double> fluxCDif_e, fluxCDif_w, fluxCDif_n, fluxCDif_s;
};

#endif //VARIABLES_DISCRETIZACION_HPP
