//
// Created by oscar on 10/09/25.
//

#include "flujo_de_masa.hpp"
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
     const Campo::Momentum                & vel,
     const Campo::Presion                 & Pstar,
     const Gradiente                      & gradP,
     const std::vector<double>            & vol,
     const Malla::Mallador::Interpolacion & inter,
     const Campo::A_coef                  & A_u,
     const Campo::A_coef                  & A_v

)
{
    coef_d.calcular(A_u, A_v);

}
