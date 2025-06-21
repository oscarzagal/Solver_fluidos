#include <utility>
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <iostream>

#include "malla_por_bloques.hpp"

namespace Malla {

  //Lista de inicializacion
  Mallador::Mallador
  (
    std::vector<int> nodos_en_x_,
    std::vector<int> nodos_en_y_,
    std::vector<std::pair<double,double>> coordenadas_en_x_,
    std::vector<std::pair<double,double>> coordenadas_en_y_,
    std::vector<std::string> frontera_norte_,
    std::vector<std::string> frontera_sur_,
    std::vector<std::string> frontera_este_,
    std::vector<std::string> frontera_oeste_
  ) :
  // std::move ayuda a mover objetos sin hacer copias innecesarias
  nodos_en_x(std::move(nodos_en_x_)),
  nodos_en_y(std::move(nodos_en_y_)),
  coordenadas_en_x(std::move(coordenadas_en_x_)),
  coordenadas_en_y(std::move(coordenadas_en_y_)),
  frontera_norte(std::move(frontera_norte_)),
  frontera_sur(std::move(frontera_sur_)),
  frontera_este(std::move(frontera_este_)),
  frontera_oeste(std::move(frontera_oeste_))
  {}

  std::vector<double> Mallador::deltas(const int n, const double h) {

    std::vector<double> delta(n);

    const double del=h/(static_cast<double>(n)-2.0);

    delta[0] = 0.0;
    delta[n-1] = 0.0;

    for (int i = 1; i < n - 1; ++i) {
      delta[i] = del;
    }

    return delta;

  }

  void Mallador::puntos
  (
      const int n, // numero de nodos por parche
      const double c_inicial,
      const double c_final,
      std::vector<double> &c_tmp
  )
  {

    std::vector<double> coord_rel(n);
    std::vector<double> delta(n);

    delta = deltas(n, c_final - c_inicial);

    coord_rel[0] = c_inicial;
    c_tmp.push_back(c_inicial);

    for (int i = 1; i < n - 1; ++i) {
      coord_rel[i] = coord_rel[i - 1] + (delta[i] + delta[i - 1]) / 2.0;
      c_tmp.push_back(coord_rel[i]);
    }

    c_tmp.push_back(c_final);
  }

  std::vector<int> Mallador::frontera_horizontal
  (
    const int lim_inf,
    const int lim_sup,
    const int nx,
    const int nodo_altura
  )
  {
    std::vector<int> parche;

    for (int i=lim_inf; i<lim_sup; ++i) {
      parche.push_back(i+nx*nodo_altura);
    }

    return parche;

  }

  std::vector<int> Mallador::frontera_vertical
  (
      const int lim_inf,
      const int lim_sup,
      const int nx,
      const int nodo_horizontal
  )
  {
    std::vector<int> parche;

    for (int j=lim_inf; j<lim_sup; ++j) {
      parche.push_back(nodo_horizontal+nx*j);
    }

    return parche;

  }

  void Mallador::cortar_nodos
  (
    const std::vector<int>& nodos,
    std::vector<int>& lim_inf_nodos,
    std::vector<int>& lim_sup_nodos
  )
  {
    /* Variable de acumulacion */
    int acum=0;

    /* Inicializacion de los nodos inferiores en cero */
    for (int i=0; i<nodos.size(); ++i) {
      lim_inf_nodos.push_back(0);
    }

    /* Asignacion de los nodos superiores al vector "lim_sup_nodos" */
    for (int j=0; j<nodos.size(); ++j) {
      for (int i=0; i<j+1; ++i) {
        acum=acum+nodos[i];
      }
      lim_sup_nodos.push_back(acum);
      acum=0;
    }

    /* Nodos cortados */
    int nodos_cortados=0;

    for (int i=0; i<nodos.size(); ++i) {
      if (i==0 || i==nodos.size()-1) {
        nodos_cortados+=1;
      }
      if (i>0 && i<nodos.size()-1) {
        nodos_cortados+=2;
      }

      lim_sup_nodos[i]=lim_sup_nodos[i]-nodos_cortados;
      if (i<nodos.size()-1) {
        lim_inf_nodos[i+1]=lim_sup_nodos[i];
      }
    }

  }

  void Mallador::comprobacion_de_tamanos
  (
    const std::vector<int> &nodos,
    const std::vector<std::pair<double, double> > &coordenadas,
    const std::vector<std::string> &nombres_frontera_1,
    const std::vector<std::string> &nombres_frontera_2,
    const std::string &direccion
  )
  {
    // Comprobacion de los tamaños de los vectores en x
    if (
      nodos.size() != coordenadas.size() ||
      coordenadas.size() != nombres_frontera_1.size() ||
      nombres_frontera_1.size() != nombres_frontera_2.size()
    ) {
      throw std::runtime_error
      (
        "La longitud de los vectores en " + direccion + " no es coincidente"
      );
    }
  }

  std::vector<double> Mallador::obtener_coordenadas_tmp_x() {

    // Limpieza del vector antes de cada llamada
    deltax.clear();

    // Coordenadas temporales
    std::vector<double> x_tmp(0.0);

    // Contenedor de nodos en x
    std::vector<int> nodos_x;

    // Coordenadas de las aristas en x
    std::vector<double> coord_inicial_x;
    std::vector<double> coord_final_x;

    // Comprobacion de los tamaños de los vectores
    comprobacion_de_tamanos
    (
      nodos_en_x,
      coordenadas_en_x,
      frontera_norte,
      frontera_sur,
      "x"
    );

    // Nodos en x
    for (const int& i : nodos_en_x) {
      nodos_x.push_back(i);
    }

    // Creacion de las aristas en x
    for (int i=0; i<coordenadas_en_x.size(); ++i) {
      coord_inicial_x.push_back(coordenadas_en_x[i].first);
      coord_final_x.push_back(coordenadas_en_x[i].second);
    }

    // Obtencion de nodos en x */
    for (int i=0; i<nodos_x.size(); ++i) {
      puntos(nodos_x[i],coord_inicial_x[i],coord_final_x[i],x_tmp);
      std::vector<double> resultado = deltas(nodos_x[i], coord_final_x[i] - coord_inicial_x[i]);
      deltax.insert(deltax.end(), resultado.begin(), resultado.end());
    }

    // Para eliminar los ceros intermedios
    deltax.erase(
      std::remove(deltax.begin() + 1, deltax.end() - 2, 0.0),deltax.end() - 2
    );

    // Nodos que deben de ser borrados en x
    for (int i=0; i<coord_final_x.size()-1; ++i) {
      x_tmp.erase
      (
        std::remove(x_tmp.begin(),x_tmp.end(),coord_final_x[i]),x_tmp.end()
      );
    }

    /* "std::remove" mueve los elementos COINCIDENTES al final del rango y
    /* retorna un iterador (un pointer) que almacena su posicion */
    /* https://en.cppreference.com/w/cpp/container/vector/begin */

    return x_tmp;

  }

  std::vector<double> Mallador::obtener_coordenadas_tmp_y() {

    // Limpieza del vector antes de cada llamada
    deltay.clear();

    // Coordenadas temporales
    std::vector<double> y_tmp(0.0);

    // Contenedor de nodos en y
    std::vector<int> nodos_y;

    // Coordenadas de las aristas en y
    std::vector<double> coord_inicial_y;
    std::vector<double> coord_final_y;

    // Comprobacion de los tamaños de los vectores
    comprobacion_de_tamanos
    (
      nodos_en_y,
      coordenadas_en_y,
      frontera_este,
      frontera_oeste,
      "y"
    );

    // Nodos en y
    for (const int& i : nodos_en_y) {
      nodos_y.push_back(i);
    }

    // Creacion de las aristas en y
    for (int i=0; i<coordenadas_en_y.size(); ++i) {
      coord_inicial_y.push_back(coordenadas_en_y[i].first);
      coord_final_y.push_back(coordenadas_en_y[i].second);
    }

    // Obtencion de nodos en y */
    for (int i=0; i<nodos_y.size(); ++i) {
      puntos(nodos_y[i],coord_inicial_y[i],coord_final_y[i],y_tmp);
      std::vector<double> resultado = deltas(nodos_y[i], coord_final_y[i] - coord_inicial_y[i]);
      deltay.insert(deltay.end(), resultado.begin(), resultado.end());
    }

    // Para eliminar los ceros intermedios
    deltay.erase(
      std::remove(deltay.begin() + 1, deltay.end() - 2, 0.0),deltay.end() - 2
    );

    // Nodos que deben de ser borrados en y
    for (int i=0; i<coord_final_y.size()-1; ++i) {
      y_tmp.erase
      (
        std::remove(y_tmp.begin(),y_tmp.end(),coord_final_y[i]),y_tmp.end()
      );
    }

    return y_tmp;

  }

  void Mallador::asignar_coord_pers_x(const std::vector<double>& x_) {
    x=x_;
  }

  void Mallador::asignar_coord_pers_y(const std::vector<double>& y_) {
    y=y_;
  }

  std::vector<double> Mallador::obtener_coord_pers_x() const {
    return x;
  }

  std::vector<double> Mallador::obtener_coord_pers_y() const {
    return y;
  }


  void Mallador::preparar_coordenadas_persistentes() {

    // Coordenadas persistentes
    std::vector<double> x;
    std::vector<double> y;

    // Coordenadas temporales
    const std::vector<double> x_tmp=this->obtener_coordenadas_tmp_x();
    const std::vector<double> y_tmp=this->obtener_coordenadas_tmp_y();

    // Numero de nodos
    const int nx=static_cast<int>(x_tmp.size());
    const int ny=static_cast<int>(y_tmp.size());

    for (int j=0; j<ny; ++j) {
      for (int i=0; i<nx; ++i) {
        x.push_back(x_tmp[i]);
        y.push_back(y_tmp[j]);
      }
    }

    // for (int j = 0; j < ny; ++j) {
    //   for (int i = 0; i < nx; ++i) {
    //     printf(
    //       "x[%d] = %f, y[%d] = %f\n", i+nx*j,
    //       x[i+nx*j], i+nx*j,
    //       y[i+nx*j]
    //     );
    //   }
    // }

    // Asignacion de coordenadas
    asignar_coord_pers_x(x);
    asignar_coord_pers_y(y);

  }

  void Mallador::asignar_parches(const std::vector<Parche> &parche, const Frontera frontera) {
    switch (frontera) {
      case Frontera::Norte:
        parches_norte=parche;
        break;
      case Frontera::Sur:
        parches_sur=parche;
        break;
      case Frontera::Este:
        parches_este=parche;
        break;
      case Frontera::Oeste:
        parches_oeste=parche;
        break;
    }
  }

  std::vector<Mallador::Parche> Mallador::obtener_parches(const Frontera frontera) const {
    switch (frontera) {
      case Frontera::Norte:
        return parches_norte;
      case Frontera::Sur:
        return parches_sur;
      case Frontera::Este:
        return parches_este;
      case Frontera::Oeste:
        return parches_oeste;
    }
    return {};
  }

  void Mallador::preparar_parches_fronteras() {

    // Lista de los nodos inferiores en x
    std::vector<int> lim_inf_nodos_x(0);

    // Lista de los nodos superiores en x
    std::vector<int> lim_sup_nodos_x(0);

    // Lista de los nodos inferiores en y
    std::vector<int> lim_inf_nodos_y(0);

    // Lista de los nodos superiores en y
    std::vector<int> lim_sup_nodos_y(0);

   // Numero de nodos por parche
    std::vector<int> nodos_x;
    std::vector<int> nodos_y;

    // Nodos en x
    for (const int& i : nodos_en_x) {
      nodos_x.push_back(i);
    }

    // Nodos en y
    for (const int& i : nodos_en_y) {
      nodos_y.push_back(i);
    }

    // Numero de nodos en x e y
    const int nx=static_cast<int>(this->obtener_coordenadas_tmp_x().size());
    const int ny=static_cast<int>(this->obtener_coordenadas_tmp_y().size());

    // La funcion "cortar_nodos" calcula los limites de las aristas en base a
    // los nodos cortados
    if (nodos_x.size()>1) {
      cortar_nodos(nodos_x,lim_inf_nodos_x,lim_sup_nodos_x);
    } else {
      lim_inf_nodos_x.push_back(0);
      lim_sup_nodos_x.push_back(nx);
    }

    if (nodos_y.size()>1) {
      cortar_nodos(nodos_y,lim_inf_nodos_y,lim_sup_nodos_y);
    } else {
      lim_inf_nodos_y.push_back(0);
      lim_sup_nodos_y.push_back(ny);
    }

    // Declaracion de los parches para la frontera norte
    std::vector<Parche> parches_norte(nodos_x.size());

    // Declaracion de los parches para la frontera sur
    std::vector<Parche> parches_sur(nodos_x.size());

    // Declaracion de los parches para la frontera este
    std::vector<Parche> parches_este(nodos_y.size());

    // Declaracion de los parches para la frontera oeste
    std::vector<Parche> parches_oeste(nodos_y.size());

    // Asignacion de los nodos frontera a los parches correspondientes
    // Al primero le corresponde la posicion 0 de lim_inf_nodos_x y
    // lim_sup_nodos_x, al segundo la posicion 1 y asi
    // NOTE: Frontera norte y sur
    for (int i=0; i<nodos_x.size(); ++i) {
      parches_norte[i].obtener_nodos_del_parche=frontera_horizontal
      (
        lim_inf_nodos_x[i],lim_sup_nodos_x[i],nx,ny-1
      );
      parches_sur[i].obtener_nodos_del_parche=frontera_horizontal
      (
        lim_inf_nodos_x[i],lim_sup_nodos_x[i],nx,0
      );
      parches_norte[i].obtener_nombre=frontera_norte[i];
      parches_sur[i].obtener_nombre=frontera_sur[i];
    }

    // NOTE: Frontera este y oeste
    for (int i=0; i<nodos_y.size(); ++i) {
      parches_este[i].obtener_nodos_del_parche=frontera_vertical
      (
        lim_inf_nodos_y[i],lim_sup_nodos_y[i],nx,nx-1
      );
      parches_oeste[i].obtener_nodos_del_parche=frontera_vertical
      (
        lim_inf_nodos_y[i],lim_sup_nodos_y[i],nx,0
      );
      parches_este[i].obtener_nombre=frontera_este[i];
      parches_oeste[i].obtener_nombre=frontera_oeste[i];
    }

    // Asignacion de los parches
    asignar_parches(parches_norte,Frontera::Norte);
    asignar_parches(parches_sur,Frontera::Sur);
    asignar_parches(parches_este,Frontera::Este);
    asignar_parches(parches_oeste,Frontera::Oeste);

  }


} // Fin namespace Mallador