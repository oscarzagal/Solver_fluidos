//
// Created by oscar on 12/06/25.
//

#ifndef ESQUEMAS_DE_DISCRETIZACION_HPP
#define ESQUEMAS_DE_DISCRETIZACION_HPP

#include "malla_por_bloques.hpp"
#include "ecuaciones_gobernantes.hpp"

namespace Esquemas_discretizacion {

    struct grad_explicito {

        // Gradiente de presion explicito multiplicado por el volumen
        std::vector<double> Pstar_x_vol, Pstar_y_vol;

    };

    void laplaciano
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
        const std::vector<double> &Pstar,
        const Malla::Mallador &malla
    );

} // Fin namespace Esquemas_discretizacion


#endif //ESQUEMAS_DE_DISCRETIZACION_HPP