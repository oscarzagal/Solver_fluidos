//
// Created by oscar on 22/06/25.
//

#include "Campo.hpp"

namespace Campo {

// Implementación de los constructores

velFace::velFace(const int nx, const int ny)
    : uFace_x(nx * ny, 0.0), vFace_y(nx * ny, 0.0)
{}

MDotStar::MDotStar(const int nx, const int ny)
    : mDotStar_x(nx * ny, 0.0), mDotStar_x_old(nx * ny, 0.0),
      mDotStar_y(nx * ny, 0.0), mDotStar_y_old(nx * ny, 0.0)
{}

Momentum::Momentum(const int nx, const int ny)
    : u_star(nx * ny, 0.0), u_old(nx * ny, 0.0),
      v_star(nx * ny, 0.0), v_old(nx * ny, 0.0),
      mdotstar(nx, ny),
      A_u{std::vector<double>(nx * ny, 0.0), std::vector<double>(nx * ny, 0.0),
          std::vector<double>(nx * ny, 0.0), std::vector<double>(nx * ny, 0.0),
          std::vector<double>(nx * ny, 0.0), std::vector<double>(nx * ny, 0.0)},
      A_v{std::vector<double>(nx * ny, 0.0), std::vector<double>(nx * ny, 0.0),
          std::vector<double>(nx * ny, 0.0), std::vector<double>(nx * ny, 0.0),
          std::vector<double>(nx * ny, 0.0), std::vector<double>(nx * ny, 0.0)}
{}

Presion::Presion(const int nx, const int ny)
    : P_star(nx * ny, 0.0), P_old(nx * ny, 0.0),
      A_p{std::vector<double>(nx * ny, 0.0), std::vector<double>(nx * ny, 0.0),
          std::vector<double>(nx * ny, 0.0), std::vector<double>(nx * ny, 0.0),
          std::vector<double>(nx * ny, 0.0), std::vector<double>(nx * ny, 0.0)}
{}


} // Fin namespace Campos
