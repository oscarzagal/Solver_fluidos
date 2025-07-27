//
// Created by oscar on 13/06/25.
//

#include<iostream>
#include<fstream>

#include "escritura.hpp"

void escribir
(
  const std::string& nombre_archivo,
  const std::string& nombre_variable,
  const std::vector<double>& x,
  const std::vector<double>& y,
  const int& nx,
  const int& ny,
  const std::vector<double>& phi
)
{

  std::fstream myFile(nombre_archivo, std::ios::out);

  myFile << "title=" << nombre_variable << std::endl;
  myFile << "variables=x,y," << nombre_variable << std::endl;
  myFile << "zone i=" << nx << "," << "j=" << ny << std::endl;
  myFile << "" << std::endl;

  for (int j=0; j<ny; ++j) {
    for (int i=0; i<nx; ++i) {
      myFile << x[i+nx*j] << " " << y[i+nx*j] << " " << phi[i+nx*j] << "\n";
    }
  }

  std::cout << "Escribiendo a " << nombre_archivo << "\n";

  myFile.close();

}
