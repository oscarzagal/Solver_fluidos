#pragma once

#include <vector>
#include <string>
#include <utility>

#include "malla_por_bloques.hpp"

/* Alias para ese tipo tan feo */
typedef std::vector<std::pair<double,double>> par_de;

typedef std::vector<int> cantidad_de;

typedef std::vector<std::string> asignacion_de;

typedef std::vector<Malla::Mallador::Parche> almacenar;

/* Nodos por cada parche en la direccion x */
inline cantidad_de nodos_en_x = {5};

/* Nodos por cada parche en la direccion y */
inline cantidad_de nodos_en_y = {5};

/* TODO: buscar la manera de tirar una excepcion ante la configuracion
 erronea en las coordenadas, por ejemplo:

    const par_de coordenadas_en_x = {
        {0.0,1.0},
        {0.0,2.0}
    };

    La segunda coordenada arranca en cero, lo cual no esta bien, es necesario
    comprobar que el "first" del segundo elemento sea igual a "second"
    del primer elemento
*/

/* Par de coordenas en x (first: coordenada inicial, second: coordenada final) */
inline par_de coordenadas_en_x = {
    {0.0, 1.0}
};

/* Par de coordenas en y (first: coordenada inicial, second: coordenada final) */
inline par_de coordenadas_en_y = {
    {0.0, 1.0}
};

/* Nombres de los parches horizontales */
inline asignacion_de nombres_frontera_norte = {
    {"norte"}
};

inline asignacion_de nombres_frontera_sur = {
    {"sur"}
};

/* Nombres de los parches verticales */
inline asignacion_de nombres_frontera_este = {
    {"este"}
};

inline asignacion_de nombres_frontera_oeste = {
    {"oeste"}
};
