#include <iostream>
#include <stdexcept>
#include <vector>

constexpr int nx = 155;
constexpr int ny = 155;

namespace Estructura {

template <typename T>
class Vector2D {
public:
  Vector2D(int nx_, int ny_, T val_inicial_) :
    nx(nx_), ny(ny_), val_inicial(val_inicial_), matriz(nx_*ny_,val_inicial_) {}

  // Sobrecarga del operador () para escritura y lectura
  T& operator()(int i, int j) {
    if (i >= nx || j >= ny) {
      throw std::out_of_range("Indice de la matriz fuera de rango");
    }
    return matriz[i+nx*j];
  }

  Vector2D<T>& operator+=(const Vector2D& otra) {
    if (nx != otra.nx || ny != otra.ny) {
      throw std::out_of_range("Los vectores son de diferente tamaño");
    }
    for (int j=0; j<ny; ++j) {
      for (int i=0; i<nx; ++i) {
        matriz[i+nx*j] += otra.matriz[i+nx*j];
      }
    }
    return *this;
  }

  // Es necesario que sea una copia para no sobreescribir a "esta"
  friend Vector2D<T> operator+(Vector2D<T> esta, const Vector2D<T>& otra) {
    esta += otra;
    return esta;
  }

  int nx;
  int ny;
  T val_inicial;
  std::vector<T> matriz;

};

}

int main() {

  Estructura::Vector2D<double> vel_u(nx,ny,3.0);
  Estructura::Vector2D<double> vel_v(nx,ny,12.0);
  Estructura::Vector2D<double> vel_w(nx,ny,10.0);

  Estructura::Vector2D<double> vel(nx,ny,0.0);

  vel = vel_u + vel_v + vel_w;

  std::cout << vel(0,0) << "\n";
  std::cout << vel_u(0,0) << "\n";
  std::cout << vel_v(0,0) << "\n";
  std::cout << vel_w(0,0) << "\n";

  return 0;
}
