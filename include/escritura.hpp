//
// Created by oscar on 13/06/25.
//

#ifndef ESCRITURA_HPP
#define ESCRITURA_HPP

#include <string>
#include <vector>

void escribir
(
  const std::string& nombre_archivo,
  const std::string& nombre_variable,
  const std::vector<double>& x,
  const std::vector<double>& y,
  const int& nx,
  const int& ny,
  const std::vector<double>& phi
);


#endif //ESCRITURA_HPP
