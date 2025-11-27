#include <algorithm>
#include <cstdio>
#include <iostream>
#include <iterator>
#include <ostream>
#include <utility>
#include <vector>
#include <fstream>

std::vector<double> deltas(const int nn, const double h) {

    std::vector<double> deltas(nn);

    const double del=h/(static_cast<double>(nn)-2.0);

    deltas[0] = 0.0;
    deltas[nn - 1] = 0.0;


    for (int i = 1; i < nn - 1; ++i) {
      deltas[i] = del;
    }


    return deltas;
}

// Modifica el estado de "centroides"
void centroides
(
    const int nn,
    const double centroide_inicial,
    const double centroide_final,
    std::vector<double> &centroides
)
{

    std::vector<double> delta = deltas(nn, centroide_final - centroide_inicial);

    centroides[0] = centroide_inicial;
    centroides[nn - 1] = centroide_final;

    for (int i = 1; i < nn - 1; ++i) {
      centroides[i] = centroides[i - 1] + (delta[i] + delta[i - 1]) / 2.0;
    }


}

void escribir
(
  const std::string& nombre_archivo,
  const std::string& nombre_variable,
  const std::vector<double> x,
  const std::vector<double> y,
  const std::vector<std::pair<int, int>> nodos,
  const std::vector<double> phi
)
{

  std::fstream myFile(nombre_archivo, std::ios::out);


  myFile << "TITLE=" << nombre_variable << std::endl;
  myFile << "variables=x,y," << nombre_variable << std::endl;

  int numero_de_elemento = 0;
  for (int i = 0 ; i < static_cast<int>(nodos.size()) ; ++i) {

    const int nx = nodos[i].first;
    const int ny = nodos[i].second;

    myFile << "zone t=\"zona" << i + 1 << "\", i=" << nx << "," << "j=" << ny << " datapacking=point" << "\n";


    for (int j = 0 ; j < ny ; ++j) {
      for (int i = 0 ; i < nx ; ++i) {
        myFile << x[numero_de_elemento] << " " << y[numero_de_elemento] << " " << phi[numero_de_elemento] << "\n";
        numero_de_elemento++;
      }
    }


  }

  std::cout << "Escribiendo a " << nombre_archivo << "\n";

  myFile.close();

}


int main () {

    // Deben de ser coincidentes los nodos que comparten caras

    const int uno_nx = 5;
    const int uno_ny = 5;

    const int dos_nx = 5;
    const int dos_ny = 5;

    const int tres_nx = 5;
    const int tres_ny = 5;

    std::vector<std::pair<int, int>> nodos = {{uno_nx, uno_ny}, {dos_nx, dos_ny}, {tres_nx, tres_ny}};

    std::vector<double> uno_centroides_x(uno_nx);
    std::vector<double> uno_centroides_y(uno_ny);

    std::vector<double> dos_centroides_x(dos_nx);
    std::vector<double> dos_centroides_y(dos_ny);

    std::vector<double> tres_centroides_x(tres_nx);
    std::vector<double> tres_centroides_y(tres_ny);

    // Coordenadas globales
    std::vector<double> x;
    std::vector<double> y;

    // Variable muda
    std::vector<double> phi((uno_nx * uno_ny) + (dos_nx * dos_ny) + (tres_nx * tres_ny), 0.0);

    // Generacion del primer bloque
    centroides(uno_nx, 0.0, 1.0, uno_centroides_x);
    centroides(uno_ny, 0.0, 1.0, uno_centroides_y);

    // Generacion del segundo bloque
    centroides(dos_nx, 1.0, 2.0, dos_centroides_x);
    centroides(dos_ny, 0.0, 1.0, dos_centroides_y);

    // Generacion del tercer bloque
    centroides(tres_nx, 0.0, 1.0, tres_centroides_x);
    centroides(tres_ny, 1.0, 2.0, tres_centroides_y);

    for (int j = 0 ; j < uno_ny ; ++j) {
      for (int i = 0 ; i < uno_nx ; ++i) {

            x.push_back(uno_centroides_x[i]);
            y.push_back(uno_centroides_y[j]);

      }
    }

    std::cout << "x.size() = " << x.size() << "\n";

    for (int j = 0 ; j < dos_ny ; ++j) {
      for (int i = 0 ; i < dos_nx ; ++i) {

            x.push_back(dos_centroides_x[i]);
            y.push_back(dos_centroides_y[j]);

      }
    }

    std::cout << "x.size() = " << x.size() << "\n";

    for (int j = 0 ; j < tres_ny ; ++j) {
      for (int i = 0 ; i < tres_nx ; ++i) {

            x.push_back(tres_centroides_x[i]);
            y.push_back(tres_centroides_y[j]);

      }
    }

    std::cout << "x.size() = " << x.size() << "\n";

    // Vector que almace las coordenadas por pares
    std::vector<std::pair<double, double>> coordenadas(x.size());

    for (int i = 0 ; i < static_cast<int>(x.size()) ; ++i) {
        coordenadas[i].first = x[i];
        coordenadas[i].second = y[i];
    }

    // TODO: utilizar 2 "std::vector<std::pair<double, double>>"
    // NOTE: A "coordenadas_tmp" se le aplica "sort" y luego "unique" para
    // obtener los elementos repetidos. Luego de eso se tiene que formar otro
    // "std::vector<std::pair<double, double>>" con los elementos duplicados,
    // para despues usarlo con "remove" y eliminar los elementos coincidentes de
    // "coordenadas" uno por uno.

    // Coordenadas temporales
    std::vector<std::pair<double, double>> coordenadas_tmp = coordenadas;

    // Ordenamiento
    std::sort(coordenadas_tmp.begin(), coordenadas_tmp.end());


    // "std::unique" para mover los duplicados al final
    auto ultimo = std::unique(coordenadas_tmp.begin(), coordenadas_tmp.end());

    // Borrar los duplicados
    coordenadas_tmp.erase(ultimo, coordenadas_tmp.end());

    // Vector de diferencia
    std::vector<std::pair<double, double>> diferencia;


    // Encontrar la diferencia
    std::set_difference(coordenadas_tmp.begin(), coordenadas_tmp.end(), coordenadas.begin(), coordenadas.end(), std::back_inserter(diferencia));

    for (const auto dif : diferencia) {

        std::cout << dif.first << " " << dif.second  << "\n";

    }

    for (const auto& dif : diferencia) {

        auto iguales = std::remove(coordenadas.begin(), coordenadas.end(), dif);

        coordenadas.erase(iguales, coordenadas.end());

    }


    // for (int i = 0 ; i < static_cast<int>(coordenadas.size()) ; ++i) {
    //     printf("x[%d] = %f y[%d] = %f\n", i, coordenadas[i].first, i, coordenadas[i].second);
    // }

    x.clear();
    y.clear();

    for (int i = 0 ; i < static_cast<int>(coordenadas.size()) ; ++i) {
        x.push_back(coordenadas[i].first);
        y.push_back(coordenadas[i].second);
    }

    std::cout << "nodos.size() = " << nodos.size() << "\n";

    escribir("malla.dat", "phi", x, y, nodos, phi);




    return 0;
}
