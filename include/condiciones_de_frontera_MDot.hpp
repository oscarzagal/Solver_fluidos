//
// Created by oscar on 28/08/25.
//

#ifndef CONDICIONES_DE_FRONTERA_MDOT_HPP
#define CONDICIONES_DE_FRONTERA_MDOT_HPP

#include <array>
#include <vector>
#include <string>

#include "config_CF.hpp"
#include "malla_por_bloques.hpp"


/*-----------------------------------------------------------------------------
                 Struct para almacenar parches Flujo de masa
-----------------------------------------------------------------------------*/

struct Parches_Flujo_de_Masa {

    // Vector que va a almacenar una copia de los nodos obtenidos de la clase
    // Malla::Mallador de los archivos "malla_por_bloques.*"
    std::vector<int> obtener_nodos_del_parche;
    std::vector<double> delta; // Ya sea en "x" o "y"
    std::string que_delta_fue_asignada;
    std::string obtener_nombre;
    std::string tipo_de_CF;
    std::string frontera_fisica;
    double vecUnitNormal; // Vector normal unitario
    int desfase;
    int nx, ny; // Nodos en "x" e "y"

    /* Funciones miembro */

    // Modifica el estado de "obtener_nodos_del_parche"
    void cortar_nodos_esquina();

    // Modifica el estado de "vecUnitNormal" y "frontera_fisica". Retorna un
    // vector normal unitario acorde a un sistema de coordenadas cartesiano
    // para la malla.
    void calcular_vector_normal_unitario();

    // Modifica el estado de "tipo_de_CF". Como los tipos de CF deben de ser
    // iguales para ambas velocidades se escoge, por convencion, los configurados
    // para "u" en "include/config_CF.hpp"
    std::string añadir_tipo_de_CF
    (
        std::array<CF_Dirichlet, limite_num_parches> g_dirichlet,
        std::array<CF_Zero_Neumann, limite_num_parches> g_zero_neumann
    );

    // Pre calcula el desfase para asignar las CF de flujo de masa en base a los
    // valores de las velocidades de los centros de masa. Modifica a "desfase"
    void calcular_desfase();

    // Modifica el estado de "delta"
    void asignar_deltas(const Malla::Mallador& malla);

    // Constructor
    Parches_Flujo_de_Masa(int, int);

};

/*-----------------------------------------------------------------------------
                 Fin Struct para almacenar parches Flujo de masa
-----------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------
                            Especializaciones
-----------------------------------------------------------------------------*/

// Template specialization para las condiciones de frontera de flujo de masa
struct Dirichlet_MDot {
    const Parches_Flujo_de_Masa& parche;
    const std::vector<double>& vel; // Velocidad en el centro de masa
    std::vector<double>& mdotstar; // Objeto a modificar
};

struct Zero_Neumann_MDot {
    const Parches_Flujo_de_Masa& parche;
    const std::vector<double>& vel; // Velocidad en el centro de masa
    std::vector<double>& mdotstar; // Objeto a modificar
};

// Clase generica (no se usa)
template<typename T>
class CF_MDot {
public:
    void aplicar() {}
};

// Especializacion para Dirichlet_MDot
template<>
class CF_MDot<Dirichlet_MDot> {
public:

    CF_MDot
    (
        const Parches_Flujo_de_Masa& parche_,
        const std::vector<double>& vel_,
        std::vector<double>& mdotstar_
    ) :
        dirichlet{parche_, vel_, mdotstar_} {}

    void aplicar() {

        const auto& nodos_parche = dirichlet.parche.obtener_nodos_del_parche;
        const auto& desfase      = dirichlet.parche.desfase;

        int index = 0;
        for (const int nodo : nodos_parche) {

            auto& mdotstar          = dirichlet.mdotstar[nodo+desfase];
            const auto& vel         = dirichlet.vel[nodo];
            const auto& delta       = dirichlet.parche.delta[index];
            const auto& vecNormUnit = dirichlet.parche.vecUnitNormal;

            mdotstar = vel * delta * vecNormUnit;

            ++index;

        }

    }

    Dirichlet_MDot dirichlet;

};

template<>
class CF_MDot<Zero_Neumann_MDot> {
public:
    CF_MDot
    (
        const Parches_Flujo_de_Masa& parche_,
        const std::vector<double>& vel_,
        std::vector<double>& mdotstar_
    ) :
        zero_neumann{parche_, vel_, mdotstar_}
    {}

    void aplicar() {

        const auto& nodos_parche = zero_neumann.parche.obtener_nodos_del_parche;
        const auto& desfase      = zero_neumann.parche.desfase;

        int index = 0;
        for (const int nodo : nodos_parche) {

            auto& mdotstar          = zero_neumann.mdotstar[nodo+desfase];
            const auto& vel         = zero_neumann.vel[nodo];
            const auto& delta       = zero_neumann.parche.delta[index];
            const auto& vecNormUnit = zero_neumann.parche.vecUnitNormal;

            mdotstar = vel * delta * vecNormUnit;

            ++index;

        }
    }

    Zero_Neumann_MDot zero_neumann;
};



/*-----------------------------------------------------------------------------
                            Fin Especializaciones
-----------------------------------------------------------------------------*/



/*-----------------------------------------------------------------------------
                        Funciones de construccion de CF
-----------------------------------------------------------------------------*/

// Modifica el estado de "lista_Dirichlet_x," "lista_Zero_Neumann_x",
// "lista_Dirichlet_y" y "lista_Zero_Neumann_y"
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

// Funcion de asignacion que va colocada dentro de "construir_CF_flujo_de_masa".
// Modifica el estado de "lista_Dirichlet" y "lista_Zero_Neumann".
void asignar_condiciones_de_frontera_MDot
(
    const std::vector<Parches_Flujo_de_Masa> & parches,
    const std::vector<double>                & vel_star,
    std::vector<double>                      & mDotStar,
    std::vector<CF_MDot<Dirichlet_MDot>>     & lista_Dirichlet,
    std::vector<CF_MDot<Zero_Neumann_MDot>>  & lista_Zero_Neumann
);


/*-----------------------------------------------------------------------------
                      Fin Funciones de construccion de CF
-----------------------------------------------------------------------------*/

#endif //CONDICIONES_DE_FRONTERA_MDOT_HPP
