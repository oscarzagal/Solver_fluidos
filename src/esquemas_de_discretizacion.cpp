//
// Created by oscar on 26/07/25.
//

#include "esquemas_de_discretizacion.hpp"
#include "config_control.hpp"
#include "utilidades.hpp"
#include <algorithm>
#include <vector>

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
            // const std::vector<double>& x = malla.x_tmp;
            // const std::vector<double>& y = malla.y_tmp;

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
                    fluxFDif_e[idx(i, j, nx)] = - gamma * deltay[j] / (x[i + 1] - x[i]);
                    fluxFDif_w[idx(i, j, nx)] = - gamma * deltay[j] / (x[i] - x[i - 1]);
                    fluxFDif_n[idx(i, j, nx)] = - gamma * deltax[i] / (y[j + 1] - y[j]);
                    fluxFDif_s[idx(i, j, nx)] = - gamma * deltax[i] / (y[j] - y[j - 1]);

                    // Coeficientes difusivos en las caras para el coeficiente agrupado a_C
                    fluxCDif_e[idx(i, j, nx)] = - fluxFDif_e[idx(i, j, nx)];
                    fluxCDif_w[idx(i, j, nx)] = - fluxFDif_w[idx(i, j, nx)];
                    fluxCDif_n[idx(i, j, nx)] = - fluxFDif_n[idx(i, j, nx)];
                    fluxCDif_s[idx(i, j, nx)] = - fluxFDif_s[idx(i, j, nx)];
                }
            }

        } // Fin funcion laplaciano_lineal

        void laplaciano_lineal_presion
        (
            const int nx,
            const int ny,
            const Coeficiente_d                  & coef_d,
            fluxes_difusivos                     & fluxes,
            const Malla::Mallador::Interpolacion & inter,
            const Malla::Mallador                & malla
        )
        {

            // Area de las caras locales
            const std::vector<double>& deltax = malla.deltax;
            const std::vector<double>& deltay = malla.deltay;

            // Coordenadas de los centroides de las celdas
            const auto& x = malla.obtener_coord_pers_x();
            const auto& y = malla.obtener_coord_pers_y();

            // Flux difusivo
            auto& fluxFDif_e = fluxes.fluxFDif_e;
            auto& fluxFDif_w = fluxes.fluxFDif_w;
            auto& fluxFDif_n = fluxes.fluxFDif_n;
            auto& fluxFDif_s = fluxes.fluxFDif_s;
            auto& fluxCDif_e = fluxes.fluxCDif_e;
            auto& fluxCDif_w = fluxes.fluxCDif_w;
            auto& fluxCDif_n = fluxes.fluxCDif_n;
            auto& fluxCDif_s = fluxes.fluxCDif_s;

            for (int j = 1 ; j < ny - 1 ; ++j) {

               // Area para las caras "e" y "w" (variable local)
               const double S_ew = deltay[j];

              for (int i = 1 ; i < nx - 1 ; ++i) {

                  /*-----------------------------------------------------------------------------
                                                Variables locales
                  -----------------------------------------------------------------------------*/

                  // Area para las caras "n" y "s"
                  const double S_ns = deltax[i];

                  const int Centro = i + nx * j;
                  const int Este   = (i + 1) + nx * j;
                  const int Oeste  = (i - 1) + nx * j;
                  const int Norte  = i + nx * (j + 1);
                  const int Sur    = i + nx * (j - 1);

                  // Coeficente "d" en los centroides de las celdas
                  const double dC_u = coef_d.dC_u[Centro];
                  const double dE_u = coef_d.dE_u[Centro];
                  const double dW_u = coef_d.dW_u[Centro];
                  const double dC_v = coef_d.dC_v[Centro];
                  const double dN_v = coef_d.dN_v[Centro];
                  const double dS_v = coef_d.dS_v[Centro];

                  // Factores de interpolacion
                  const double ge = inter.ge[Centro];
                  const double gw = inter.gw[Centro];
                  const double gn = inter.gn[Centro];
                  const double gs = inter.gs[Centro];

                  // Coeficiente "d" inteprolado
                  const double de_i = interpolar(dC_u, dE_u, ge);
                  const double dw_i = interpolar(dC_u, dW_u, gw);
                  const double dn_i = interpolar(dC_v, dN_v, gn);
                  const double ds_i = interpolar(dC_v, dS_v, gs);

                  // "δx_{CE}", "δx_{CW}", "δy_{CN}" y "δy_{CS}"
                  const double delta_x_CE = x[Este] - x[Centro];
                  const double delta_x_CW = x[Centro] - x[Oeste];
                  const double delta_y_CN = y[Norte] - y[Centro];
                  const double delta_y_CS = y[Centro] - y[Sur];

                  /*-----------------------------------------------------------------------------
                                             Fin Variables locales
                  -----------------------------------------------------------------------------*/

                  // Coeficientes difusivos en las caras para el coeficiente agrupado a_F
                  fluxFDif_e[Centro] = - de_i * S_ew / delta_x_CE;
                  fluxFDif_w[Centro] = - dw_i * S_ew / delta_x_CW;
                  fluxFDif_n[Centro] = - dn_i * S_ns / delta_y_CN;
                  fluxFDif_s[Centro] = - ds_i * S_ns / delta_y_CS;

                  // Coeficientes difusivos en las caras para el coeficiente agrupado a_C
                  fluxCDif_e[Centro] = fluxFDif_e[Centro];
                  fluxCDif_w[Centro] = fluxFDif_w[Centro];
                  fluxCDif_n[Centro] = fluxFDif_n[Centro];
                  fluxCDif_s[Centro] = fluxFDif_s[Centro];

              }
            }


        } // Fin funcion laplaciano_lineal_presion

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
            const auto &mDotStar_x = mstar.mDotStar_x;
            const auto &mDotStar_y = mstar.mDotStar_y;

            for (int j = 1 ; j < ny - 1 ; ++j) {
                for (int i = 1 ; i< nx - 1 ; ++i) {

                    // Renombrar de nuevo por legibilidad
                    const auto& mDotStar_e = mDotStar_x[idx(i+1,j,nx)];
                    const auto& mDotStar_w = mDotStar_x[idx(i,j,nx)];
                    const auto& mDotStar_n = mDotStar_y[idx(i,j+1,nx)];
                    const auto& mDotStar_s = mDotStar_y[idx(i,j,nx)];

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

                for (int i = 1 ; i < nx - 1 ; ++i) {

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

                    // NOTE: recordar que los gradientes de presion van multiplicados por
                    // el inverso de la densidad en la ecuacion de momentum
                    gPstar_x_vol[Centro] = (P_e - P_w) * Deltay;
                    gPstar_y_vol[Centro] = (P_n - P_s) * Deltax;
                }
            }

        }

    } // Fin namespace Explicita


    void construccion_coeficiente_b_momemtum
    (
         const int nx,
         const int ny,
         Campo::A_coef             & A,
         const std::vector<double> & vol,
         const std::vector<double> & grad_vol,
         const std::vector<double> & vel
    )
    {

        auto& b        = A.b;
        const auto& ac = A.ac;

        const double inv_rho = 1.0 / rho;

        for (int j = 1 ; j < ny - 1 ; ++j) {
            for (int i = 1 ; i < nx - 1 ; ++i) {

                const int Centro = i + nx * j;

                const auto Grad_vol = grad_vol[Centro];
                const auto Ac       = ac[Centro];
                const auto Vel      = vel[Centro];
                const auto Vol      = vol[Centro];

                /*                                 discretizacion implicita      +   pseudo transitorio */
                b[Centro] = - inv_rho * Grad_vol + (1.0 - lambda_Vel) * Ac * Vel + Vel * Vol / delta_t;

            }
        }

    }


    void construccion_matriz_A_momentum
    (
         const int nx,
         const int ny,
         const fluxes_difusivos    & fluxes_dif,
         const fluxes_convectivos  & fluxes_conv,
         const std::vector<double> & vol,
         const std::vector<double> & vel_u,
         const std::vector<double> & vel_v,
         Campo::A_coef             & A_u,
         Campo::A_coef             & A_v,
         Gradiente                 & grad
    )
    {

        auto& ae = A_u.ae;
        auto& aw = A_u.aw;
        auto& an = A_u.an;
        auto& as = A_u.as;
        auto& ac = A_u.ac;

        for (int j = 1 ; j < ny - 1 ; ++j) {
            for (int i = 1 ; i < nx - 1 ; ++i) {

                const int Centro = i + nx * j;

                /*-----------------------------------------------------------------------------
                  Coeficientes vecinos
                  -----------------------------------------------------------------------------*/

                const auto fluxFDif_e = fluxes_dif.fluxFDif_e[Centro];
                const auto fluxFDif_w = fluxes_dif.fluxFDif_w[Centro];
                const auto fluxFDif_n = fluxes_dif.fluxFDif_n[Centro];
                const auto fluxFDif_s = fluxes_dif.fluxFDif_s[Centro];

                const auto fluxFConv_e = fluxes_conv.fluxFConv_e[Centro];
                const auto fluxFConv_w = fluxes_conv.fluxFConv_w[Centro];
                const auto fluxFConv_n = fluxes_conv.fluxFConv_n[Centro];
                const auto fluxFConv_s = fluxes_conv.fluxFConv_s[Centro];

                /*-----------------------------------------------------------------------------
                  Fin Coeficientes vecinos
                  -----------------------------------------------------------------------------*/



                /*-----------------------------------------------------------------------------
                  Coeficientes centrales
                  -----------------------------------------------------------------------------*/

                const auto fluxCDif_e = fluxes_dif.fluxCDif_e[Centro];
                const auto fluxCDif_w = fluxes_dif.fluxCDif_w[Centro];
                const auto fluxCDif_n = fluxes_dif.fluxCDif_n[Centro];
                const auto fluxCDif_s = fluxes_dif.fluxCDif_s[Centro];

                const auto fluxCConv_e = fluxes_conv.fluxCConv_e[Centro];
                const auto fluxCConv_w = fluxes_conv.fluxCConv_w[Centro];
                const auto fluxCConv_n = fluxes_conv.fluxCConv_n[Centro];
                const auto fluxCConv_s = fluxes_conv.fluxCConv_s[Centro];

                const auto sumFluxCDif  = fluxCDif_e + fluxCDif_w + fluxCDif_n + fluxCDif_s;
                const auto sumFluxCConv = fluxCConv_e + fluxCConv_w + fluxCConv_n + fluxCConv_s;

                /*-----------------------------------------------------------------------------
                  Fin Coeficientes centrales
                  -----------------------------------------------------------------------------*/

                ae[Centro] = fluxFDif_e + fluxFConv_e;
                aw[Centro] = fluxFDif_w + fluxFConv_w;
                an[Centro] = fluxFDif_n + fluxFConv_n;
                as[Centro] = fluxFDif_s + fluxFConv_s;
                ac[Centro] = (sumFluxCDif + sumFluxCConv + vol[Centro] / delta_t) / lambda_Vel;

            }
        }

        // Asignacion para los coeficientes de "v", recordar que son los mismos,
        // el unico que cambia es "b".
        A_v.ae = ae;
        A_v.aw = aw;
        A_v.an = an;
        A_v.as = as;
        A_v.ac = ac;

        Discretizacion::construccion_coeficiente_b_momemtum(nx, ny, A_u, vol, grad.grad_x_vol, vel_u);
        Discretizacion::construccion_coeficiente_b_momemtum(nx, ny, A_v, vol, grad.grad_y_vol, vel_v);

    } // Fin funcion construccion_matriz_A_momentum


    void construccion_matriz_A_presion
    (
        const int nx,
        const int ny,
        const fluxes_difusivos & fluxes_difusivos,
        const Campo::MDotStar  & mdotstar,
        Campo::A_coef          & A_p
    )
    {

        auto& ae = A_p.ae;
        auto& aw = A_p.aw;
        auto& an = A_p.an;
        auto& as = A_p.as;
        auto& ac = A_p.ac;
        auto& b  = A_p.b;

        for (int j = 1 ; j < ny - 1 ; ++j) {
            for (int i = 1 ; i < nx - 1 ; ++i) {

                const int Centro = i + nx * j;
                const int Este   = (i + 1) + nx * j;
                const int Norte  = i + nx * (j + 1);

                /*-----------------------------------------------------------------------------
                  Coeficientes vecinos
                  -----------------------------------------------------------------------------*/

                const double fluxFDif_e = fluxes_difusivos.fluxFDif_e[Centro];
                const double fluxFDif_w = fluxes_difusivos.fluxFDif_w[Centro];
                const double fluxFDif_n = fluxes_difusivos.fluxFDif_n[Centro];
                const double fluxFDif_s = fluxes_difusivos.fluxFDif_s[Centro];

                /*-----------------------------------------------------------------------------
                  Fin Coeficientes vecinos
                  -----------------------------------------------------------------------------*/

                /*-----------------------------------------------------------------------------
                  Coeficientes centrales
                  -----------------------------------------------------------------------------*/

                const double fluxCDif_e = fluxes_difusivos.fluxCDif_e[Centro];
                const double fluxCDif_w = fluxes_difusivos.fluxCDif_w[Centro];
                const double fluxCDif_n = fluxes_difusivos.fluxCDif_n[Centro];
                const double fluxCDif_s = fluxes_difusivos.fluxCDif_s[Centro];

                const double sumFluxCDif = fluxCDif_e + fluxCDif_w + fluxCDif_n + fluxCDif_s;

                /*-----------------------------------------------------------------------------
                  Fin Coeficientes centrales
                  -----------------------------------------------------------------------------*/

                /*-----------------------------------------------------------------------------
                  Flujo de masa
                  -----------------------------------------------------------------------------*/

                const double mDotStar_e = mdotstar.mDotStar_x[Este];
                const double mDotStar_w = mdotstar.mDotStar_x[Centro];
                const double mDotStar_n = mdotstar.mDotStar_y[Norte];
                const double mDotStar_s = mdotstar.mDotStar_y[Centro];

                const double sumMDotStar = mDotStar_e + mDotStar_w + mDotStar_n + mDotStar_s;


                /*-----------------------------------------------------------------------------
                  Fin Flujo de masa
                  -----------------------------------------------------------------------------*/

                ae[Centro] = fluxFDif_e;
                aw[Centro] = fluxFDif_w;
                an[Centro] = fluxFDif_n;
                as[Centro] = fluxFDif_s;
                ac[Centro] = - sumFluxCDif;
                b[Centro]  = - sumMDotStar;

            }
        }

    } // Fin funcion construccion_matriz_A_presion


} // Fin namespace Esquemas_discretizacion
