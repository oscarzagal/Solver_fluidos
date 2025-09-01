//
// Created by oscar on 28/08/25.
//

#ifndef CONDICIONES_DE_FRONTERA_MDOT_HPP
#define CONDICIONES_DE_FRONTERA_MDOT_HPP

#include "malla_por_bloques.hpp"
#include <vector>
#include <string>

/*-----------------------------------------------------------------------------
                            Especializaciones
-----------------------------------------------------------------------------*/

// Template specialization para las condiciones de frontera de flujo de masa
struct Dirichlet_MDot {
    const std::vector<int>& nodos_del_parche;
    const double valor;
    std::vector<double>& mdotstar;
};

struct Zero_Neumann_MDot {
    std::vector<double>& mdotstar;
    const std::vector<int>& nodos_del_parche;
    const std::string& frontera_fisica;
    const int nx;
};

// Clase generica (no se usa)
template<typename T>
class CF_MDot {
public:
    void aplicar() {}
};

// Especializacion para Dirichlet_MDot
// TODO: Implementar la aplicacion de la condicion de frontera
template<>
class CF_MDot<Dirichlet_MDot> {
public:

    CF_MDot<Dirichlet_MDot>
    (
        // TODO: se van a requerir los campos de velocidad, entonces es necesario quitar a "valor_"
        const std::vector<int>& nodos_,
        const double valor_,
        std::vector<double>& mdotstar_
    ) :
        dirichlet{nodos_, valor_, mdotstar_}
    {}

    void aplicar() {
        // NOTE: se tiene para consulta
        // for (const auto& nodo : dirichlet.nodos_del_parche) {
        //     dirichlet.mdotstar[nodo] = dirichlet.valor;
        // }
    }

    Dirichlet_MDot dirichlet;

};

template<>
class CF_MDot<Zero_Neumann_MDot> {
public:
};

/*-----------------------------------------------------------------------------
                            Fin Especializaciones
-----------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------
                 Struct para almacenar parches Flujo de masa
-----------------------------------------------------------------------------*/

struct Parches_Flujo_de_Masa {

    // Vector que va a almacenar una copia de los nodos obtenidos de la clase
    // Malla::Mallador de los archivos "malla_por_bloques.*"
    std::vector<int> obtener_nodos_del_parche;
    std::string obtener_nombre;
    double vecUnitNormal; // Vector normal unitario
    int nx, ny; // Nodos en "x" e "y"

    /* Funciones miembro */

    // Modifica el estado de "obtener_nodos_del_parche"
    void cortar_nodos_esquina();

    // Modifica el estado de "vecUnitNormal". Retorna un vector normal unitario
    // acorde a un sistema de coordenadas cartesiano para la malla.
    void calcular_vector_normal_unitario();

    // Constructor
    Parches_Flujo_de_Masa(int, int);

};

/*-----------------------------------------------------------------------------
                 Fin Struct para almacenar parches Flujo de masa
-----------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------
                        Funciones de construccion de CF
-----------------------------------------------------------------------------*/

// TODO: cambiar malla por los parches especiales para flujo de masa
void construir_CF_flujo_de_masa
(
    const Malla::Mallador                   & malla,
    const std::vector<double>               & u_star,
    const std::vector<double>               & v_star,
    std::vector<double>                     & mDotStar_x,
    std::vector<double>                     & mDotStar_y,
    std::vector<CF_MDot<Dirichlet_MDot>>    & lista_Dirichlet,
    std::vector<CF_MDot<Zero_Neumann_MDot>> & lista_Zero_Neumann
);

// TODO: Implementar la funcion "asignar_condiciones_de_frontera_MDot" similar a la de "asignar_condiciones_de_frontera" en "condiciones_de_frontera.*"
// Usar "emplace_back" en lugar de "push_back" para evitar temporales innecesarios. Más información en "tests/template_specialization.cpp"


/*-----------------------------------------------------------------------------
                      Fin Funciones de construccion de CF
-----------------------------------------------------------------------------*/

#endif //CONDICIONES_DE_FRONTERA_MDOT_HPP
