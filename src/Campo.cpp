//
// Created by oscar on 22/06/25.
//

#include "Campo.hpp"

namespace Campo {

// Implementación de los constructores

velFace::velFace(const int nx, const int ny, double campo_inicial)
    : uFace_x(nx * ny, campo_inicial), vFace_y(nx * ny, campo_inicial),
      uFace_x_n(nx * ny, campo_inicial), vFace_y_n(nx * ny, campo_inicial)
{}

MDotStar::MDotStar(const int nx, const int ny, double campo_inicial)
    : mDotStar_x(nx * ny, campo_inicial), mDotStar_x_old(nx * ny, campo_inicial),
      mDotStar_y(nx * ny, campo_inicial), mDotStar_y_old(nx * ny, campo_inicial)
{}

Momentum::Momentum(const int nx, const int ny, double campo_inicial_u, double campo_inicial_v)
    : u_star(nx * ny, campo_inicial_u), u_old(nx * ny, campo_inicial_u),
      v_star(nx * ny, campo_inicial_v), v_old(nx * ny, campo_inicial_v),
      A_u{std::vector<double>(nx * ny), std::vector<double>(nx * ny),
          std::vector<double>(nx * ny), std::vector<double>(nx * ny),
          std::vector<double>(nx * ny), std::vector<double>(nx * ny)},
      A_v{std::vector<double>(nx * ny), std::vector<double>(nx * ny),
          std::vector<double>(nx * ny), std::vector<double>(nx * ny),
          std::vector<double>(nx * ny), std::vector<double>(nx * ny)}
{}

Presion::Presion(const int nx, const int ny, double campo_inicial)
    : P_star(nx * ny, campo_inicial), P_old(nx * ny, campo_inicial),
      A_p{std::vector<double>(nx * ny), std::vector<double>(nx * ny),
          std::vector<double>(nx * ny), std::vector<double>(nx * ny),
          std::vector<double>(nx * ny), std::vector<double>(nx * ny)}
{}


} // Fin namespace Campos
