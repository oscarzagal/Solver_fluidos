//
// Created by oscar on 26/07/25.
//

#include "esquemas_de_discretizacion.hpp"
#include "utilidades.hpp"

namespace Esquemas_discretizacion {

    void laplaciano
    (
        const int &nx,
        const int &ny,
        const double &k,
        Ecuaciones_gobernantes::A_coef &A,
        Malla::Mallador &malla
    ) {
        /* Obtencion de las deltas */
        // Las deltas dependen del tamaño de la malla y de los nodos
        std::vector<double> deltax = malla.deltax;
        std::vector<double> deltay = malla.deltay;

        /* Obtencion de las cordenadas */
        std::vector<double> x = malla.obtener_coordenadas_tmp_x();
        std::vector<double> y = malla.obtener_coordenadas_tmp_y();

        A.ae.resize(nx * ny);
        A.aw.resize(nx * ny);
        A.an.resize(nx * ny);
        A.as.resize(nx * ny);
        A.ap.resize(nx * ny);
        A.b.resize(nx * ny);

        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                A.ae[i + nx * j] = k * deltay[j] / (x[i + 1] - x[i]);
                A.aw[i + nx * j] = k * deltay[j] / (x[i] - x[i - 1]);
                A.an[i + nx * j] = k * deltax[i] / (y[j + 1] - y[j]);
                A.as[i + nx * j] = k * deltax[i] / (y[j] - y[j - 1]);
                A.ap[i + nx * j] = A.ae[i + nx * j] + A.aw[i + nx * j] + A.an[i + nx * j] + A.as[i + nx * j];
                A.b[i + nx * j] = 0.0;
            }
        }
    }

    void gradiente_explicito
    (
        const int nx,
        const int ny,
        const Malla::Mallador::Interpolacion &inter,
        grad_explicito &grad_explicito,
        const std::vector<double> &Pstar,
        const Malla::Mallador &malla
    )
    {
        grad_explicito.Pstar_x_vol.resize(nx * ny);
        grad_explicito.Pstar_y_vol.resize(nx * ny);

        std::vector<double> deltax = malla.deltax;
        std::vector<double> deltay = malla.deltay;

        auto &gPstar_x_vol = grad_explicito.Pstar_x_vol;
        auto &gPstar_y_vol = grad_explicito.Pstar_y_vol;

        const auto &ge = inter.ge;
        const auto &gw = inter.gw;
        const auto &gn = inter.gn;
        const auto &gs = inter.gs;

        for (int j=1; j<ny-1; ++j) {
            for (int i=1; i<nx-1; ++i) {
                // TODO: terminar
                gPstar_x_vol[idx(i,j,nx)]=(ge[idx(i,j,nx)]*Pstar[idx(i+1,j,nx)]
                    +(1.0-ge[idx(i,j,nx)])*Pstar[idx(i,j,nx)]-gw[idx(i,j,nx)]
                    *Pstar[idx(i-1,j,nx)]-(1.0-gw[idx(i,j,nx)])*Pstar[idx(i,j,nx)])
                    *deltay[j];
                gPstar_y_vol[idx(i,j,nx)]=
            }
        }

    }


} // Fin namespace Esquemas_discretizacion