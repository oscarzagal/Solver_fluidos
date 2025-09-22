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

// NOTE: asignarles el valor incial de 1.0 a los coeficientes agrupados es necesario para evitar
// divisiones sobre cero en `flujo_de_masa.cpp::Coeficiente_d::calcular`. Pero la justificacion
// más importante está relacionada con las condiciones de frontera de Dirichlet y Zero Neumann
// para la velocidad y la presion: para ambos tipos el coeficiente central será 1.0. La
// explicacion del porque tiene que ver con la ecuacion generadora:
//
//  a_{C}\phi_{C}+\sum_{F \in \{E,W,N,S,T,B\}}a_{F}\phi_{F}=b_{C}
//
// Se aprecia que, al hacer a_{C} = 1.0 y también al valor correspondiente en la frontera a_{F},
// (por ejemplo, este), mientras que todos los demas 0.0 se tiene la condicion de frontera
// Zero Neumann. Si a_{C} = 1.0 y b_{C} tiene un valor especificado por el usuario, mientras que los
// demas son cero, se tiene la condicion de frontera de Dirichlet.
//
// Debido a como está construido el codigo en cuestion de las condiciones de frontera para
// variables escalares, solo es necesario que a_{C} = 1.0 debido a la razon expuesta arriba.
// Los otros coeficientes no son necesarios ya que en `condiciones_de_frontera.cpp` se obtienen
// los valores de frontera sin usar la ecuacion generadora, solo el correspondiente campo \phi

Momentum::Momentum(const int nx, const int ny, double campo_inicial_u, double campo_inicial_v)
    : u_star(nx * ny, campo_inicial_u), u_old(nx * ny, campo_inicial_u),
      v_star(nx * ny, campo_inicial_v), v_old(nx * ny, campo_inicial_v),
      A_u{std::vector<double>(nx * ny, 1.0), std::vector<double>(nx * ny, 0.0),
          std::vector<double>(nx * ny, 0.0), std::vector<double>(nx * ny, 0.0),
          std::vector<double>(nx * ny, 0.0), std::vector<double>(nx * ny, 0.0)},
      A_v{std::vector<double>(nx * ny, 1.0), std::vector<double>(nx * ny, 0.0),
          std::vector<double>(nx * ny, 0.0), std::vector<double>(nx * ny, 0.0),
          std::vector<double>(nx * ny, 0.0), std::vector<double>(nx * ny, 0.0)}
{}

Presion::Presion(const int nx, const int ny, double campo_inicial)
    : P_star(nx * ny, campo_inicial), P_old(nx * ny, campo_inicial),
      A_p{std::vector<double>(nx * ny, 1.0), std::vector<double>(nx * ny, 0.0),
          std::vector<double>(nx * ny, 0.0), std::vector<double>(nx * ny, 0.0),
          std::vector<double>(nx * ny, 0.0), std::vector<double>(nx * ny, 0.0)}
{}


} // Fin namespace Campos
