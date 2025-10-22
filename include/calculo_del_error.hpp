//
// Created by oscar on 13/06/25.
//

#ifndef CALCULO_DEL_ERROR_HPP
#define CALCULO_DEL_ERROR_HPP

#include <vector>

double calcular_error_mayor
(
  const int& nx,
  const int& ny,
  const std::vector<double>& phi,
  const std::vector<double>& phi_old
);

#endif //CALCULO_DEL_ERROR_HPP
