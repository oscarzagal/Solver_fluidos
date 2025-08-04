//
// Created by oscar on 26/07/25.
//

#include "esquemas_de_discretizacion.hpp"
#include "utilidades.hpp"
#include <algorithm>

namespace Esquemas_discretizacion {

void laplaciano_lineal
(
    const int &nx,
    const int &ny,
    const double &k,
    Ecuaciones_gobernantes::A_coef &A,
    Malla::Mallador &malla
) {
    // Obtencion de las deltas
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
    Ecuaciones_gobernantes::Momentum::mstar &mstar
)
{

    // Cambiar de tamaño
    fluxes.fluxFConv_e.resize(nx*ny);
    fluxes.fluxFConv_w.resize(nx*ny);
    fluxes.fluxFConv_n.resize(nx*ny);
    fluxes.fluxFConv_s.resize(nx*ny);
    fluxes.fluxCConv_e.resize(nx*ny);
    fluxes.fluxCConv_w.resize(nx*ny);
    fluxes.fluxCConv_n.resize(nx*ny);
    fluxes.fluxCConv_s.resize(nx*ny);
    mstar.me_star.resize(nx*ny);
    mstar.mw_star.resize(nx*ny);
    mstar.mn_star.resize(nx*ny);
    mstar.ms_star.resize(nx*ny);

    // Renombrar por legibilidad
    auto &fluxFConv_e = fluxes.fluxFConv_e;
    auto &fluxFConv_w = fluxes.fluxFConv_w;
    auto &fluxFConv_n = fluxes.fluxFConv_n;
    auto &fluxFConv_s = fluxes.fluxFConv_s;
    auto &fluxCConv_e = fluxes.fluxCConv_e;
    auto &fluxCConv_w = fluxes.fluxCConv_w;
    auto &fluxCConv_n = fluxes.fluxCConv_n;
    auto &fluxCConv_s = fluxes.fluxCConv_s;
    const auto &me_star = mstar.me_star;
    const auto &mw_star = mstar.mw_star;
    const auto &mn_star = mstar.mn_star;
    const auto &ms_star = mstar.ms_star;

    for (int j=1; j<ny-1; ++j) {
        for (int i=1; i<nx-1; ++i) {
            // Flux para las celdas vecinas
            fluxFConv_e[idx(i,j,nx)] = -std::max(0.0,-me_star[idx(i,j,nx)]);
            fluxFConv_w[idx(i,j,nx)] = -std::max(0.0,-mw_star[idx(i,j,nx)]);
            fluxFConv_n[idx(i,j,nx)] = -std::max(0.0,-mn_star[idx(i,j,nx)]);
            fluxFConv_s[idx(i,j,nx)] = -std::max(0.0,-ms_star[idx(i,j,nx)]);

            // Flux para la celda central
            fluxCConv_e[idx(i,j,nx)] = std::max(0.0,me_star[idx(i,j,nx)]);
            fluxCConv_w[idx(i,j,nx)] = std::max(0.0,mw_star[idx(i,j,nx)]);
            fluxCConv_n[idx(i,j,nx)] = std::max(0.0,mn_star[idx(i,j,nx)]);
            fluxCConv_s[idx(i,j,nx)] = std::max(0.0,ms_star[idx(i,j,nx)]);
        }
    }

}


} // Fin namespace Esquemas_discretizacion
