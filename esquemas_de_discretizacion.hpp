//
// Created by oscar on 12/06/25.
//

#ifndef ESQUEMAS_DE_DISCRETIZACION_HPP
#define ESQUEMAS_DE_DISCRETIZACION_HPP

#include "malla_por_bloques.hpp"

namespace Esquemas_discretizacion {

    template<typename T>
    void laplaciano
    (
        const int &nx,
        const int &ny,
        const double &k,
        T &c,
        Malla::Mallador& malla
    ) {
        /* Obtencion de las deltas */
        // Las deltas dependen del tamaño de la malla y de los nodos
        std::vector<double> deltax = malla.deltax;
        std::vector<double> deltay = malla.deltay;

        /* Obtencion de las cordenadas */
        std::vector<double> x = malla.obtener_coordenadas_tmp_x();
        std::vector<double> y = malla.obtener_coordenadas_tmp_y();

        c.ae.resize(nx * ny);
        c.aw.resize(nx * ny);
        c.an.resize(nx * ny);
        c.as.resize(nx * ny);
        c.ap.resize(nx * ny);
        c.b.resize(nx * ny);

        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                c.ae[i + nx * j] = k * deltay[j] / (x[i + 1] - x[i]);
                c.aw[i + nx * j] = k * deltay[j] / (x[i] - x[i - 1]);
                c.an[i + nx * j] = k * deltax[i] / (y[j + 1] - y[j]);
                c.as[i + nx * j] = k * deltax[i] / (y[j] - y[j - 1]);
                c.ap[i + nx * j] = c.ae[i + nx * j] + c.aw[i + nx * j] + c.an[i + nx * j] + c.as[i + nx * j];
                c.b[i + nx * j] = 0.0;
            }
        }
    }
} // Fin namespace Esquemas_discretizacion


#endif //ESQUEMAS_DE_DISCRETIZACION_HPP
