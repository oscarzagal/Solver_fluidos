//
// Created by oscar on 26/07/25.
//

#include "esquemas_de_discretizacion.hpp"
#include "utilidades.hpp"
#include <algorithm>

namespace Esquemas_discretizacion {

void laplaciano_lineal
(
    const int nx,
    const int ny,
    const double gamma,
    fluxes_difusivos &fluxes,
    Malla::Mallador &malla
) {
    // Obtencion de las deltas
    // Las deltas dependen del tamaño de la malla y de los nodos
    std::vector<double> deltax = malla.deltax;
    std::vector<double> deltay = malla.deltay;

    // Obtencion de las cordenadas
    std::vector<double> x = malla.obtener_coordenadas_tmp_x();
    std::vector<double> y = malla.obtener_coordenadas_tmp_y();

    // Asignar tamaño
    fluxes.fluxFDif_e.resize(nx*ny,0.0);
    fluxes.fluxFDif_w.resize(nx*ny,0.0);
    fluxes.fluxFDif_n.resize(nx*ny,0.0);
    fluxes.fluxFDif_s.resize(nx*ny,0.0);
    fluxes.fluxCDif_e.resize(nx*ny,0.0);
    fluxes.fluxCDif_w.resize(nx*ny,0.0);
    fluxes.fluxCDif_n.resize(nx*ny,0.0);
    fluxes.fluxCDif_s.resize(nx*ny,0.0);

    // Renombrar por legibilidad
    auto &fluxFDif_e = fluxes.fluxFDif_e;
    auto &fluxFDif_w = fluxes.fluxFDif_w;
    auto &fluxFDif_n = fluxes.fluxFDif_n;
    auto &fluxFDif_s = fluxes.fluxFDif_s;
    auto &fluxCDif_e = fluxes.fluxCDif_e;
    auto &fluxCDif_w = fluxes.fluxCDif_w;
    auto &fluxCDif_n = fluxes.fluxCDif_n;
    auto &fluxCDif_s = fluxes.fluxCDif_s;


    for (int j = 1; j < ny - 1; j++) {
        for (int i = 1; i < nx - 1; i++) {
            // Coeficientes difusivos en las caras para el coeficiente agrupado a_F
            fluxFDif_e[idx(i,j,nx)] = gamma * deltay[j] / (x[i + 1] - x[i]);
            fluxFDif_w[idx(i,j,nx)] = gamma * deltay[j] / (x[i] - x[i - 1]);
            fluxFDif_n[idx(i,j,nx)] = gamma * deltax[i] / (y[j + 1] - y[j]);
            fluxFDif_s[idx(i,j,nx)] = gamma * deltax[i] / (y[j] - y[j - 1]);
        }
    }

    // Coeficientes difusivos en las caras para el coeficiente agrupado a_C
    fluxCDif_e = fluxFDif_e;
    fluxCDif_w = fluxFDif_w;
    fluxCDif_n = fluxFDif_n;
    fluxCDif_s = fluxFDif_s;
}

void gradiente_explicito
(
    const int nx,
    const int ny,
    const Malla::Mallador::Interpolacion &inter,
    Grad_explicito &grad_explicito,
    const std::vector<double> &P,
    const Malla::Mallador &malla
)
{
    grad_explicito.grad_x_vol.resize(nx * ny);
    grad_explicito.grad_y_vol.resize(nx * ny);

    std::vector<double> deltax = malla.deltax;
    std::vector<double> deltay = malla.deltay;

    auto &gPstar_x_vol = grad_explicito.grad_x_vol;
    auto &gPstar_y_vol = grad_explicito.grad_y_vol;

    const auto &ge = inter.ge;
    const auto &gw = inter.gw;
    const auto &gn = inter.gn;
    const auto &gs = inter.gs;

    for (int j=1; j<ny-1; ++j) {
        for (int i=1; i<nx-1; ++i) {
            gPstar_x_vol[idx(i,j,nx)]=(interpolar(P[idx(i,j,nx)],P[idx(i+1,j,nx)],ge[idx(i,j,nx)])-
                interpolar(P[idx(i,j,nx)],P[idx(i-1,j,nx)],gw[idx(i,j,nx)]))*deltay[j];
            gPstar_y_vol[idx(i,j,nx)]=(interpolar(P[idx(i,j,nx)],P[idx(i,j+1,nx)],gn[idx(i,j,nx)])-
                interpolar(P[idx(i,j,nx)],P[idx(i,j-1,nx)],gs[idx(i,j,nx)]))*deltax[i];
        }
    }

}

void divergencia_upwind
(
    const int nx,
    const int ny,
    fluxes_convectivos &fluxes,
    Ecuaciones_gobernantes::Momentum::Mstar &mstar
)
{

    // Asignar tamaño
    fluxes.fluxFConv_e.resize(nx*ny);
    fluxes.fluxFConv_w.resize(nx*ny);
    fluxes.fluxFConv_n.resize(nx*ny);
    fluxes.fluxFConv_s.resize(nx*ny);
    fluxes.fluxCConv_e.resize(nx*ny);
    fluxes.fluxCConv_w.resize(nx*ny);
    fluxes.fluxCConv_n.resize(nx*ny);
    fluxes.fluxCConv_s.resize(nx*ny);

    // Renombrar por legibilidad
    auto &fluxFConv_e = fluxes.fluxFConv_e;
    auto &fluxFConv_w = fluxes.fluxFConv_w;
    auto &fluxFConv_n = fluxes.fluxFConv_n;
    auto &fluxFConv_s = fluxes.fluxFConv_s;
    auto &fluxCConv_e = fluxes.fluxCConv_e;
    auto &fluxCConv_w = fluxes.fluxCConv_w;
    auto &fluxCConv_n = fluxes.fluxCConv_n;
    auto &fluxCConv_s = fluxes.fluxCConv_s;
    const auto &mstar_x = mstar.mstar_x;
    const auto &mstar_y = mstar.mstar_y;

    for (int j=1; j<ny-1; ++j) {
        for (int i=1; i<nx-1; ++i) {
            // Flux para las celdas vecinas
            fluxFConv_e[idx(i,j,nx)] = -std::max(0.0,-mstar_x[idx(i+1,j,nx)]);
            fluxFConv_w[idx(i,j,nx)] = -std::max(0.0,-mstar_x[idx(i,j,nx)]);
            fluxFConv_n[idx(i,j,nx)] = -std::max(0.0,-mstar_y[idx(i,j+1,nx)]);
            fluxFConv_s[idx(i,j,nx)] = -std::max(0.0,-mstar_y[idx(i,j,nx)]);

            // Flux para la celda central
            fluxCConv_e[idx(i,j,nx)] = std::max(0.0,mstar_x[idx(i+1,j,nx)]);
            fluxCConv_w[idx(i,j,nx)] = std::max(0.0,mstar_x[idx(i,j,nx)]);
            fluxCConv_n[idx(i,j,nx)] = std::max(0.0,mstar_y[idx(i,j+1,nx)]);
            fluxCConv_s[idx(i,j,nx)] = std::max(0.0,mstar_y[idx(i,j,nx)]);
        }
    }

}


} // Fin namespace Esquemas_discretizacion
