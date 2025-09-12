//
// Created by oscar on 26/07/25.
//

#include "esquemas_de_discretizacion.hpp"
#include "utilidades.hpp"
#include <algorithm>

namespace Discretizacion {

    namespace Implicita {

        void laplaciano_lineal
        (
             const int nx,
             const int ny,
             const double gamma,
             fluxes_difusivos &fluxes,
             const Malla::Mallador &malla
        )
        {
            // Obtencion de las deltas
            // Las deltas dependen del tamaño de la malla y de los nodos
            std::vector<double> deltax = malla.deltax;
            std::vector<double> deltay = malla.deltay;

            // Obtencion de las cordenadas
            std::vector<double> x = malla.retornar_coordanas_tmp(Malla::Nodos::nx);
            std::vector<double> y = malla.retornar_coordanas_tmp(Malla::Nodos::ny);

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
                    fluxFDif_e[idx(i,j,nx)] = - gamma * deltay[j] / (x[i + 1] - x[i]);
                    fluxFDif_w[idx(i,j,nx)] = - gamma * deltay[j] / (x[i] - x[i - 1]);
                    fluxFDif_n[idx(i,j,nx)] = - gamma * deltax[i] / (y[j + 1] - y[j]);
                    fluxFDif_s[idx(i,j,nx)] = - gamma * deltax[i] / (y[j] - y[j - 1]);

                    // Coeficientes difusivos en las caras para el coeficiente agrupado a_C
                    fluxCDif_e[idx(i, j, nx)] = - fluxFDif_e[idx(i, j, nx)];
                    fluxCDif_w[idx(i, j, nx)] = - fluxFDif_w[idx(i, j, nx)];
                    fluxCDif_n[idx(i, j, nx)] = - fluxFDif_n[idx(i, j, nx)];
                    fluxCDif_s[idx(i, j, nx)] = - fluxFDif_s[idx(i, j, nx)];
                }
            }

        } // Fin funcion laplaciano_lineal

        void divergencia_upwind
        (
             const int nx,
             const int ny,
             fluxes_convectivos &fluxes,
             Campo::MDotStar &mstar
        )
        {

            // Renombrar por legibilidad
            auto &fluxFConv_e = fluxes.fluxFConv_e;
            auto &fluxFConv_w = fluxes.fluxFConv_w;
            auto &fluxFConv_n = fluxes.fluxFConv_n;
            auto &fluxFConv_s = fluxes.fluxFConv_s;
            auto &fluxCConv_e = fluxes.fluxCConv_e;
            auto &fluxCConv_w = fluxes.fluxCConv_w;
            auto &fluxCConv_n = fluxes.fluxCConv_n;
            auto &fluxCConv_s = fluxes.fluxCConv_s;
            const auto &mDotstar_x = mstar.mDotStar_x;
            const auto &mDotstar_y = mstar.mDotStar_y;

            for (int j = 1 ; j < ny - 1 ; ++j) {
                for (int i = 1 ; i< nx - 1 ; ++i) {

                    // Renombrar de nuevo por legibilidad
                    const auto& mDotStar_e = mDotstar_x[idx(i+1,j,nx)];
                    const auto& mDotStar_w = mDotstar_x[idx(i,j,nx)];
                    const auto& mDotStar_n = mDotstar_y[idx(i,j+1,nx)];
                    const auto& mDotStar_s = mDotstar_y[idx(i,j,nx)];

                    // Flux para las celdas vecinas
                    fluxFConv_e[idx(i,j,nx)] = - std::max(0.0, - mDotStar_e);
                    fluxFConv_w[idx(i,j,nx)] = - std::max(0.0, - mDotStar_w);
                    fluxFConv_n[idx(i,j,nx)] = - std::max(0.0, - mDotStar_n);
                    fluxFConv_s[idx(i,j,nx)] = - std::max(0.0, - mDotStar_s);

                    // Flux para la celda central
                    fluxCConv_e[idx(i,j,nx)] = std::max(0.0, mDotStar_e);
                    fluxCConv_w[idx(i,j,nx)] = std::max(0.0, mDotStar_w);
                    fluxCConv_n[idx(i,j,nx)] = std::max(0.0, mDotStar_n);
                    fluxCConv_s[idx(i,j,nx)] = std::max(0.0, mDotStar_s);
                }
            }

        } // Fin funcion divergencia_upwind

    } // Fin namespace Implicita


    namespace Explicita {

        void gradiente
        (
             const int nx,
             const int ny,
             const Malla::Mallador::Interpolacion &inter,
             Gradiente &grad_explicito,
             const std::vector<double> &P,
             const Malla::Mallador &malla
        )
        {

            std::vector<double> deltax = malla.deltax;
            std::vector<double> deltay = malla.deltay;

            auto &gPstar_x_vol = grad_explicito.grad_x_vol;
            auto &gPstar_y_vol = grad_explicito.grad_y_vol;


            for (int j=1; j<ny-1; ++j) {

                // Solo se calcula una vez por cada paso de "j"
                const auto Deltay = deltay[j];

                for (int i=1; i<nx-1; ++i) {

                    const int Centro = i + nx * j;
                    const int Este   = (i + 1) + nx * j;
                    const int Oeste  = (i - 1) + nx * j;
                    const int Norte  = i + nx * (j + 1);
                    const int Sur    = i + nx * (j - 1);

                    const auto ge = inter.ge[Centro];
                    const auto gw = inter.gw[Centro];
                    const auto gn = inter.gn[Centro];
                    const auto gs = inter.gs[Centro];

                    const auto P_C = P[Centro];
                    const auto P_E = P[Este];
                    const auto P_W = P[Oeste];
                    const auto P_N = P[Norte];
                    const auto P_S = P[Sur];

                    const auto Deltax = deltax[i];

                    // Se devuelve por valor porque es una copia muy barata
                    const auto P_e = interpolar(P_C, P_E, ge);
                    const auto P_w = interpolar(P_C, P_W, gw);
                    const auto P_n = interpolar(P_C, P_N, gn);
                    const auto P_s = interpolar(P_C, P_S, gs);

                    gPstar_x_vol[Centro] = (P_e - P_w) * Deltay;
                    gPstar_y_vol[Centro] = (P_n - P_s) * Deltax;
                }
            }

        }

    } // Fin namespace Explicita





} // Fin namespace Esquemas_discretizacion
