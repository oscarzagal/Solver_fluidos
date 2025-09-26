//
// Created by oscar on 10/09/25.
//

#include "flujo_de_masa.hpp"
#include "config_control.hpp"
#include "utilidades.hpp"
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

          // FIXME: se me ocurrio que se pueden ahorrar muchas variables
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
    const auto& y = malla.obtener_coord_pers_y();

    // Recordar que las deltas representan el area de las caras de los elementos
    // computacionales
    const auto& deltay = malla.deltay;
    const auto& deltax = malla.deltax;

    // Obtencion de los coeficientes "d" centrales para obtener "d_interp_x" y
    // "d_interp_y"
    coef_d.calcular(A_u, A_v);

    /*-----------------------------------------------------------------------------
                        Flujo de masa en la direccion "x"
    -----------------------------------------------------------------------------*/

    for (int j = 1 ; j < ny - 1 ; ++j) {

        // Area de la cara local "e"
        const double S_e = deltay[j];

        for (int i = 1 ; i < nx - 2 ; ++i) {

            /*-----------------------------------------------------------------------------
                                  Definicion de variables locales
            -----------------------------------------------------------------------------*/

            const int Centro = i + nx * j;
            const int Este   = (i + 1) + nx * j;

            // Factor de interpolacion
            const double gx = inter.ge[Centro];

            // Velocidades en los centroides de los elementos locales "C" y "E"
            const double uC_star = vel.u_star[Centro];
            const double uE_star = vel.u_star[Este];

            // Coeficientes "d" centrales
            const double dC_u = coef_d.dC_u[Centro];
            const double dE_u = coef_d.dE_u[Centro];

            // Gradiente de presion
            const double gradPC_x = gradP.grad_x_vol[Centro] / vol[Centro];
            const double gradPE_x = gradP.grad_x_vol[Este] / vol[Este];

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

            // Flujo de masa
            double& mdotstar_e = mdotstar.mDotStar_x[Este];


            /*-----------------------------------------------------------------------------
                                Fin Definicion de variables locales
            -----------------------------------------------------------------------------*/



            /*-----------------------------------------------------------------------------
                                Calculo del flujo de masa en la cara "e"
            -----------------------------------------------------------------------------*/

            // Rhie-Chow
            ue = uface_este_interp - d_interp_x * (gradP_e - gradP_e_interp)
                // Termino adicional para evitar la dependencia en el factor de relajacion
                + (1.0 - lambda_Vel) * (ue_n - uface_este_interp_n);

            mdotstar_e = ue * S_e;

            /*-----------------------------------------------------------------------------
                                Fin Calculo del flujo de masa en la cara "e"
            -----------------------------------------------------------------------------*/

        }
    } // Fin bucle for anidado

    /*-----------------------------------------------------------------------------
                        Fin Flujo de masa en la direccion "x"
    -----------------------------------------------------------------------------*/



    /*-----------------------------------------------------------------------------
                        Flujo de masa en la direccion "y"
    -----------------------------------------------------------------------------*/


    for (int j = 1 ; j < ny - 2 ; ++j) {
        for (int i = 1 ; i < nx - 1 ; ++i) {

            /*-----------------------------------------------------------------------------
                                Definicion de variables locales
            -----------------------------------------------------------------------------*/

            // Area de la cara local "n"
            const double S_n = deltax[i];

            const int Centro = i + nx * j;
            const int Norte  = i + nx * (j + 1);

            // Factor de interpolacion
            const double gy = inter.gn[Centro];

            // Velocidades en los centroides de los elementos locales "C" y "N"
            const double vC_star = vel.v_star[Centro];
            const double vN_star = vel.v_star[Norte];

            // Coeficientes "d" centrales
            const double dC_v = coef_d.dC_v[Centro];
            const double dN_v = coef_d.dN_v[Centro];

            // Gradiente de presion
            const double gradPC_y = gradP.grad_y_vol[Centro] / vol[Centro];
            const double gradPN_y = gradP.grad_y_vol[Norte] / vol[Norte];

            // Calculo de la velocidad interpolada en la cara local "n"
            // (se requiere modificar el miembro de "velFace", por tal motivo
            // se usa una referencia).
            double& vface_norte_interp = velFace.vFace_y_interp[Norte];
            vface_norte_interp = interpolar(vC_star, vN_star, gy);

            // Velocidad interpolada de la iteracion anterior
            // (no hace falta modificar al miembro, por lo que se usa una copia)
            const double vface_norte_interp_n = velFace.vFace_y_interp_n[Norte];

            // Coeficiente "d" interpolado
            const double d_interp_y = interpolar(dC_v, dN_v, gy);

            // Gradiente de presion en la cara interpolado
            const double gradP_n_interp = interpolar(gradPC_y, gradPN_y, gy);

            // Presiones
            const double P_C = Pstar.P_star[Centro];
            const double P_N = Pstar.P_star[Norte];

            // δX_{CN}
            const double delta_y_CN = y[Norte] - y[Centro];

            // Gradiente de presion sobre la cara "n"
            const double gradP_n = (P_N - P_C) / delta_y_CN;

            // Velocidad en la cara "n"
            // (se requiere modificar el miembro de "velFace", por tal motivo
            // se usa una referencia).
            double& vn = velFace.vFace_y[Norte];

            // Velocidad en la cara "n" de la iteracion anterior
            // (no hace falta modificar al miembro, por lo que se usa una copia)
            const double vn_n = velFace.vFace_y_n[Norte];

            // Flujo de masa
            double& mdotstar_n = mdotstar.mDotStar_y[Norte];

            /*-----------------------------------------------------------------------------
                                 Fin Definicion de variables locales
            -----------------------------------------------------------------------------*/



            /*-----------------------------------------------------------------------------
                                Calculo del flujo de masa en la cara "n"
              -----------------------------------------------------------------------------*/

            // Rhie-Chow
            vn = vface_norte_interp - d_interp_y * (gradP_n - gradP_n_interp)
                // Termino adicional para evitar la dependencia en el factor de relajacion
                + (1.0 - lambda_Vel) * (vn_n - vface_norte_interp_n);

            mdotstar_n = vn * S_n;

            /*-----------------------------------------------------------------------------
                                Fin Calculo del flujo de masa en la cara "n"
              -----------------------------------------------------------------------------*/

        }
    }

    /*-----------------------------------------------------------------------------
                    Fin Flujo de masa en la direccion "y"
    -----------------------------------------------------------------------------*/

}
