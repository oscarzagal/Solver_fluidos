//
// Created by oscar on 10/09/25.
//

#include "flujo_de_masa.hpp"
#include "config_control.hpp"
#include "utilidades.hpp"
#include <asm-generic/errno.h>
#include <vector>

// Constructor
Coeficiente_d::Coeficiente_d(const int nx_, const int ny_, const std::vector<double>& vol_)
    : dC_u(nx_ * ny_, 0.0), dE_u(nx_ * ny_, 0.0), dW_u(nx_ * ny_, 0.0), dC_v(nx_ * ny_, 0.0),
      dN_v(nx_ * ny_, 0.0), dS_v(nx_ * ny_, 0.0), vol(vol_), nx(nx_), ny(ny_) {}

void Coeficiente_d::calcular(const Campo::A_coef& A_u, const Campo::A_coef& A_v) {

    const auto& ac_u = A_u.ac;
    const auto& ac_v = A_v.ac;

    for (int j = 1 ; j < ny - 1 ; ++j) {
      for (int i = 1 ; i < nx - 1 ; ++i) {

          const int Centro = i + nx * j;
          const int Este   = (i + 1) + nx * j;
          const int Oeste  = (i - 1) + nx * j;
          const int Norte  = i + nx * (j + 1);
          const int Sur    = i + nx * (j - 1);

          const double volC = vol[Centro];
          const double volE = vol[Este];
          const double volW = vol[Oeste];
          const double volN = vol[Norte];
          const double volS = vol[Sur];

          const double acC_u = ac_u[Centro];
          const double acE_u = ac_u[Este];
          const double acW_u = ac_u[Oeste];
          const double acC_v = ac_v[Centro];
          const double acN_v = ac_v[Norte];
          const double acS_v = ac_v[Sur];

          // NOTE: se me ocurrio que se pueden ahorrar muchas variables
          // si solo se calculan los coeficientes "dC_u" y "dC_v"
          dC_u[Centro] = volC / acC_u;
          dE_u[Centro] = volE / acE_u;
          dW_u[Centro] = volW / acW_u;
          dC_v[Centro] = volC / acC_v;
          dN_v[Centro] = volN / acN_v;
          dS_v[Centro] = volS / acS_v;

      }
    }
}


void calcular_flujo_de_masa
(
     int nx,
     int ny,
     Campo::velFace                       & velFace,
     Campo::MDotStar                      & mdotstar,
     Coeficiente_d                        & coef_d,
     const Malla::Mallador                & malla,
     const Campo::Momentum                & vel,
     const Campo::Presion                 & Pstar,
     const Gradiente                      & gradP,
     const std::vector<double>            & vol,
     const Malla::Mallador::Interpolacion & inter,
     const Campo::A_coef                  & A_u,
     const Campo::A_coef                  & A_v

)
{
    // Coordenadas de los centroides de las celdas
    const auto& x = malla.obtener_coord_pers_x();

    // Obtencion de los coeficientes "d" centrales para obtener "d_interp_x" y
    // "d_interp_y"
    coef_d.calcular(A_u, A_v);

    // std::vector<double> d_interp_x(nx * ny, 0.0);
    // std::vector<double> d_interp_y(nx * ny, 0.0);

    /*-----------------------------------------------------------------------------
                        Flujo de masa en la direccion "x"
    -----------------------------------------------------------------------------*/

    for (int j = 1 ; j < ny - 1 ; ++j) {
        for (int i = 1 ; i < nx - 2 ; ++i) {

            /*-----------------------------------------------------------------------------
                                  Definicion de variables locales
            -----------------------------------------------------------------------------*/

            const int Centro = i + nx * j;
            const int Este   = (i + 1) + nx * j;

            // Factor de interpolacion
            const auto gx = inter.ge[Centro];

            // Velocidades en los centroides de los elementos locales "C" y "E"
            const auto uC_star = vel.u_star[Centro];
            const auto uE_star = vel.u_star[Este];

            // Coeficientes "d" centrales
            const auto dC_u = coef_d.dC_u[Centro];
            const auto dE_u = coef_d.dE_u[Centro];

            // Gradiente de presion por el volumen
            const double gradPC_x_times_vol = gradP.grad_x_vol[Centro];
            const double gradPE_x_times_vol = gradP.grad_x_vol[Este];

            // Volumen
            const double volC = vol[Centro];
            const double volE = vol[Este];

            // Gradiente de presion
            const double gradPC_x = gradPC_x_times_vol / volC;
            const double gradPE_x = gradPE_x_times_vol / volE;

            // Calculo de la velocidad interpolada en la cara local "e"
            // (se requiere modificar el miembro de "velFace", por tal motivo
            // se usa una referencia).
            double& uface_este_interp = velFace.uFace_x_interp[Este];
            uface_este_interp = interpolar(uC_star, uE_star, gx);

            // Velocidad interpolada de la iteracion anterior
            // (no hace falta modificar al miembro, por lo que se usa una copia)
            const double uface_este_interp_n = velFace.uFace_x_interp_n[Este];

            // Coeficiente "d" interpolado
            const double d_interp_x = interpolar(dC_u, dE_u, gx);

            // Gradiente de presion en la cara interpolado
            const double gradP_e_interp = interpolar(gradPC_x, gradPE_x, gx);

            // Presiones
            const double P_C = Pstar.P_star[Centro];
            const double P_E = Pstar.P_star[Este];

            // δX_{CE}
            const double delta_x_CE = x[Este] - x[Centro];

            // Gradiente de presion sobre la cara "e"
            const double gradP_e = (P_E - P_C) / delta_x_CE;

            // Velocidad en la cara "e"
            // (se requiere modificar el miembro de "velFace", por tal motivo
            // se usa una referencia).
            double& ue = velFace.uFace_x[Este];

            // Velocidad en la cara "e" de la iteracion anterior
            // (no hace falta modificar al miembro, por lo que se usa una copia)
            const double ue_n = velFace.uFace_x_n[Este];


            /*-----------------------------------------------------------------------------
                                Fin Definicion de variables locales
            -----------------------------------------------------------------------------*/



            /*-----------------------------------------------------------------------------
                                Calculo de la velocidad en la cara "e"
            -----------------------------------------------------------------------------*/

            ue = uface_este_interp - d_interp_x * (gradP_e - gradP_e_interp) // Rhie-Chow
                + (1.0 - lambda_Vel) * (ue_n - uface_este_interp_n); // Termino adicional


            // TODO: calcular el flujo de masa

            /*-----------------------------------------------------------------------------
                            Fin Calculo de la velocidad en la cara "e"
            -----------------------------------------------------------------------------*/






        }
    }

    /*-----------------------------------------------------------------------------
                        Fin Flujo de masa en la direccion "x"
    -----------------------------------------------------------------------------*/



    /*-----------------------------------------------------------------------------
                        Flujo de masa en la direccion "y"
    -----------------------------------------------------------------------------*/


    for (int j = 1 ; j < ny - 2 ; ++j) {
        for (int i = 1 ; i < nx - 1 ; ++i) {

            const int Centro = i + nx * j;
            const int Norte  = i + nx * (j + 1);


        }
    }

    /*-----------------------------------------------------------------------------
                    Fin Flujo de masa en la direccion "y"
    -----------------------------------------------------------------------------*/



}
