//
// Created by oscar on 08/10/25.
//

#include "convergencia.hpp"
#include "config_control.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iterator>
#include <stdexcept>



std::pair<double, std::string> error_mayor
(
    const int nx,
    const int ny,
    const Campo::Momentum & velU,
    const Campo::Presion  & presion,
    const Campo::MDotStar & mdotstar,
    std::vector<double>   & error_mayor_por_campo
)
{

    // Campos nuevos
    const auto& u_star     = velU.u_star;
    const auto& v_star     = velU.v_star;
    const auto& P_star     = presion.P_star;
    const auto& mdotstar_x = mdotstar.mDotStar_x;
    const auto& mdotstar_y = mdotstar.mDotStar_y;

    // Campos viejos
    const auto& u_old          = velU.u_old;
    const auto& v_old          = velU.v_old;
    const auto& P_old          = presion.P_old;
    const auto& mdotstar_x_old = mdotstar.mDotStar_x_old;
    const auto& mdotstar_y_old = mdotstar.mDotStar_y_old;

    error_por_campo(nx, ny, Residual::presion, P_star, P_old, error_mayor_por_campo);
    error_por_campo(nx, ny, Residual::u, u_star, u_old, error_mayor_por_campo);
    error_por_campo(nx, ny, Residual::v, v_star, v_old, error_mayor_por_campo);
    error_por_campo(nx, ny, Residual::mdotstar_x, mdotstar_x, mdotstar_x_old, error_mayor_por_campo);
    error_por_campo(nx, ny, Residual::mdotstar_y, mdotstar_y, mdotstar_y_old, error_mayor_por_campo);

    if (error_mayor_por_campo.empty()) throw std::runtime_error("El vector 'error_mayor_por_campo' esta vacio, panzón");

    const auto max_iterator = std::max_element(error_mayor_por_campo.begin(), error_mayor_por_campo.end());

    if (!std::isfinite(*max_iterator)) throw std::runtime_error("Residual de infinito o NAN detectado, panzón");

    const double mayor_error = *max_iterator;
    const int indice_del_error_mayor = std::distance(error_mayor_por_campo.begin(), max_iterator);

    static const std::array<std::string, NUM_CAMPOS> campos {
        "Presion", "u", "v", "mdotstar_x", "mdotstar_y"
    };

    return {mayor_error, campos[indice_del_error_mayor]};

}


void error_por_campo
(
    const int nx,
    const int ny,
    Residual residual,
    const std::vector<double> & phi_new,
    const std::vector<double> & phi_old,
    std::vector<double>       & error_mayor_por_campo
)
{

    double residual_local = 0.0;

    error_mayor_por_campo[static_cast<int>(residual)] = 0.0;
    for (int j = 0 ; j < ny ; ++j) {
      for (int i = 0 ; i < nx ; ++i) {

          const int Centro = i + nx * j;

          residual_local = std::abs(phi_new[Centro] - phi_old[Centro]);

          error_mayor_por_campo[static_cast<int>(residual)] =
              std::max(error_mayor_por_campo[static_cast<int>(residual)], residual_local);
      }
    }

}
