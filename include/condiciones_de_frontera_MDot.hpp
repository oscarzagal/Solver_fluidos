//
// Created by oscar on 28/08/25.
//

#ifndef CONDICIONES_DE_FRONTERA_MDOT_HPP
#define CONDICIONES_DE_FRONTERA_MDOT_HPP

#include <array>
#include <vector>
#include <string>

#include "config_CF.hpp"

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
    // TODO: Implementar.
};

/*-----------------------------------------------------------------------------
                            Fin Especializaciones
-----------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------
                 Struct para almacenar parches Flujo de masa
-----------------------------------------------------------------------------*/

// TODO: Agregaer una funcion miembro para obtener el tipo de parche (Dirichlet_MDot o
// Zero_Neumann_MDot)
struct Parches_Flujo_de_Masa {

    // Vector que va a almacenar una copia de los nodos obtenidos de la clase
    // Malla::Mallador de los archivos "malla_por_bloques.*"
    std::vector<int> obtener_nodos_del_parche;
    std::string obtener_nombre;
    std::string tipo_de_CF;
    double vecUnitNormal; // Vector normal unitario
    int nx, ny; // Nodos en "x" e "y"

    /* Funciones miembro */

    // Modifica el estado de "obtener_nodos_del_parche"
    void cortar_nodos_esquina();

    // Modifica el estado de "vecUnitNormal". Retorna un vector normal unitario
    // acorde a un sistema de coordenadas cartesiano para la malla.
    void calcular_vector_normal_unitario();

    // Modifica el estado de "tipo_de_CF". Como los tipos de CF deben de ser
    // iguales para ambas velocidades se escoge, por convencion, los configurados
    // para "u" en "include/config_CF.hpp"
    std::string añadir_tipo_de_CF
    (
        std::array<CF_Dirichlet, limite_num_parches> g_dirichlet,
        std::array<CF_Zero_Neumann, limite_num_parches> g_zero_neumann
    );

    // Constructor
    Parches_Flujo_de_Masa(int, int);

};

/*-----------------------------------------------------------------------------
                 Fin Struct para almacenar parches Flujo de masa
-----------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------
                        Funciones de construccion de CF
-----------------------------------------------------------------------------*/

// TODO: Implementar esta fucnion
void construir_CF_flujo_de_masa
(
    const std::vector<Parches_Flujo_de_Masa> & parches_norte_FM,
    const std::vector<Parches_Flujo_de_Masa> & parches_sur_FM,
    const std::vector<Parches_Flujo_de_Masa> & parches_este_FM,
    const std::vector<Parches_Flujo_de_Masa> & parches_oeste_FM,
    const std::vector<double>                & u_star,
    const std::vector<double>                & v_star,
    std::vector<double>                      & mDotStar_x,
    std::vector<double>                      & mDotStar_y,
    std::vector<CF_MDot<Dirichlet_MDot>>     & lista_Dirichlet_x,
    std::vector<CF_MDot<Zero_Neumann_MDot>>  & lista_Zero_Neumann_x,
    std::vector<CF_MDot<Dirichlet_MDot>>     & lista_Dirichlet_y,
    std::vector<CF_MDot<Zero_Neumann_MDot>>  & lista_Zero_Neumann_y
);

// Funcion de asignacion que va colocada dentro de "construir_CF_flujo_de_masa"
void asignar_condiciones_de_frontera_MDot
(
    const std::vector<Parches_Flujo_de_Masa> & parches,
    const std::vector<double>                & vel_star,
    std::vector<double>                      & mDotStar,
    std::vector<CF_MDot<Dirichlet_MDot>>     & lista_Dirichlet,
    std::vector<CF_MDot<Zero_Neumann_MDot>>  & lista_Zero_Neumann
);

// TODO: Implementar la funcion "asignar_condiciones_de_frontera_MDot" similar a la de "asignar_condiciones_de_frontera" en "condiciones_de_frontera.*"
// Usar "emplace_back" en lugar de "push_back" para evitar temporales innecesarios. Más información en "tests/template_specialization.cpp"


/*-----------------------------------------------------------------------------
                      Fin Funciones de construccion de CF
-----------------------------------------------------------------------------*/

#endif //CONDICIONES_DE_FRONTERA_MDOT_HPP
