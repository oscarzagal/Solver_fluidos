//
// Created by oscar on 1/06/25.
//
/* # MALLADOR POR BLOQUES DE REGIONES RECTANGULARES.
 A través del archivo "config_malla.hpp" es posible crear mallas rectangulares
 con un numero arbitrario de parches, lo cual es util para resolver problemas
 sencillos de dinamica de fluidos computacional.

 ## Limitaciones:
 - Solo produce mallas ortogonales, lo que reduce su uso a geometrias muy
   simples.
 - No es posible crear regiones compuestas por varios rectangulos. Esta es una
   caracteristica deseable que espero algun dia pueda implementar.

 ## Instrucciones de uso:

 Los nodos se asignan de la siguiente manera (suponiendo 3 parches):

      *
      |
  ny3 |
      |
      *
      |
  ny2 |
      |
      *
      |
  ny1 |
      |
      o-----*-----*-----*
        nx1   nx2   nx3

    y+
    ^
    |
    o-> x+

  Para asignar las coordenadas se requieren los dos puntos que componen una
  arista:

        arista
    *------------*
    ci           cf

    ci: coordenada inicial
    cf: coordenada inicial

  Se tendran tantas coordenadas como aristas.

  La asignación de los nombres de los parches debe de ser de izquierda a derecha
  para las fronteras horizontales y de abajo a arriba para las fronteras
  verticales

  TODO: escribir ejemplos

 */

#ifndef MALLA_POR_BLOQUES_HPP // guarda
#define MALLA_POR_BLOQUES_HPP

#include <utility>
#include <vector>
#include <string>

namespace Malla {

  enum class Frontera {
    Norte,
    Sur,
    Este,
    Oeste
  };

  class Mallador {


  public:

    // Constructor
    Mallador
    (
      std::vector<int>,
      std::vector<int>,
      std::vector<std::pair<double,double>>,
      std::vector<std::pair<double,double>>,
      std::vector<std::string>,
      std::vector<std::string>,
      std::vector<std::string>,
      std::vector<std::string>
    );

    // Obtencion de deltas
    static std::vector<double> deltas(int n, double h);

    // Obtencion de los puntos por parche
    static void puntos (
        int n,
        double c_inicial,
        double c_final,
        std::vector<double> &c_tmp
    );

    // Definicion de los nodos globales para las fronteras horizontales
    static std::vector<int> frontera_horizontal
    (
      int lim_inf,
      int lim_sup,
      int nx,
      int nodo_altura
    );

    // Definicion de los nodos globales para las fronteras verticales
    static std::vector<int> frontera_vertical
    (
      int lim_inf,
      int lim_sup,
      int nx,
      int nodo_horizontal
    );

    //Funcion que establece los limites de los parches en cuestion de nodos
    static void cortar_nodos
    (
      const std::vector<int>& nodos,
      std::vector<int>& lim_inf_nodos,
      std::vector<int>& lim_sup_nodos
    );

    // Funcion para comprobar los tamaños de los vectores
    static void comprobacion_de_tamanos
    (
      const std::vector<int>& nodos,
      const std::vector<std::pair<double,double>>& coordenadas,
      const std::vector<std::string>& nombres_frontera_1,
      const std::vector<std::string>& nombres_frontera_2,
      const std::string& direccion
    );

    // Coordenadas temporales para x
    // En estas dos funciones esa donde se llamaran a todas las funciones
    // estaticas, pudiendo entonces obtener y retornar los vecotres x e y
    [[nodiscard]] std::vector<double> obtener_coordenadas_tmp_x();

    // [[nodiscard]] sirve para enviar un warning si el valor retornado no es usado
    // https://en.cppreference.com/w/cpp/language/attributes/nodiscard

    // Coordenadas temporales para y
    [[nodiscard]] std::vector<double> obtener_coordenadas_tmp_y();

    [[nodiscard]] std::vector<double> obtener_volumenes();


    // Asignacion de coordenadas persistentes
    void preparar_coordenadas_persistentes();


    // Obtener coordenada persistente en x
    [[nodiscard]] std::vector<double> obtener_coord_pers_x() const;

    // Obtener coordenada persistente en y
    [[nodiscard]] std::vector<double> obtener_coord_pers_y() const;

    // Asignacion de parches frontera
    void preparar_parches_fronteras();

    // TODO: hacer un setter y un getter para "nx" y "ny"

    struct Parche {

      // Nodos que corresponden con el parche asociado
      std::vector<int> obtener_nodos_del_parche;

      std::string obtener_nombre;

    };

    struct Interpolacion {

      // Factores de interpolacion para las caras de las celdas del elemnto C
      std::vector<double> ge,gw,gn,gs;

    };

    [[nodiscard]] std::vector<Parche> obtener_parches(Frontera frontera) const;

    [[nodiscard]] static Interpolacion obtener_factores_de_interpolacion(Mallador& malla);

    // Vector de deltas para su posterior uso en los esquemas de discretizacion
    std::vector<double> deltax;
    std::vector<double> deltay;


  private:

    // Asignar coordenada persistente en x
    void asignar_coord_pers_x(const std::vector<double>& x_);

    // Asignar coordenada persistente en y
    void asignar_coord_pers_y(const std::vector<double>& y_);

    // Asignar los parches fronteras
    void asignar_parches(const std::vector<Parche>& parche, Frontera frontera);

    // Nodos por cada parche en la direccion x
    std::vector<int> nodos_en_x;

    // Nodos por cada parche en la direccion y
    std::vector<int> nodos_en_y;

    // Par de coordenas en x (first: coordenada inicial, second: coordenada final)
    std::vector<std::pair<double,double>> coordenadas_en_x;

    // Par de coordenas en y (first: coordenada inicial, second: coordenada final)
    std::vector<std::pair<double,double>> coordenadas_en_y;

    // Parches horizontales (nombre, tipo de condicion de frontera)
    std::vector<std::string> frontera_norte;
    std::vector<std::string> frontera_sur;

    // Nombres de los parches verticales
    std::vector<std::string> frontera_este;
    std::vector<std::string> frontera_oeste;

    // Coordenadas persistentes
    std::vector<double> x;
    std::vector<double> y;

    // Parches por frontera
    std::vector<Parche> parches_norte;
    std::vector<Parche> parches_sur;
    std::vector<Parche> parches_este;
    std::vector<Parche> parches_oeste;

  };

} // Fin namespace Mallador

#endif //MALLA_POR_BLOQUES_HPP
