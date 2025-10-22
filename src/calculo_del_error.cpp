//
// Created by oscar on 13/06/25.
//

#include "calculo_del_error.hpp"
#include <cmath>

double calcular_error_mayor
(
  const int& nx,
  const int& ny,
  const std::vector<double>& phi,
  const std::vector<double>& phi_old
)
{

    double error_mayor=0.0;

    for (int j=1; j<ny-1; ++j) {
        for (int i=1; i<nx-1; ++i) {
            error_mayor=std::max(error_mayor,fabs(phi[i+nx*j]-phi_old[i+nx*j]));
        }
    }

    return error_mayor;

}
