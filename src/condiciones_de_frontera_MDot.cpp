//
// Created by oscar on 31/08/25.
//

#include "condiciones_de_frontera_MDot.hpp"
#include <stdexcept>

// Constructor
Parches_Flujo_de_Masa::Parches_Flujo_de_Masa(const int nx_, const int ny_)
: nx(nx_), ny(ny_) {}

void Parches_Flujo_de_Masa::cortar_nodos_esquina() {

    if (obtener_nodos_del_parche.empty()) {
        throw std::runtime_error("El parche " + obtener_nombre + "esta vacío");
    }

    for (int i = 0; i < static_cast<int>(obtener_nodos_del_parche.size()); ++i) {

        // Se obtienen los indices locales "i" y "j"
        int index_i = obtener_nodos_del_parche[i] % nx;
        int index_j = obtener_nodos_del_parche[i] / nx;

        if ((index_i == 0      && index_j == 0)      ||
            (index_i == nx - 1 && index_j == ny - 1) ||
            (index_i == nx - 1 && index_j == 0)      ||
            (index_i == 0      && index_j == ny - 1)) {

            obtener_nodos_del_parche.erase(obtener_nodos_del_parche.begin() + i);

        }
    }

}

void Parches_Flujo_de_Masa::calcular_vector_normal_unitario() {

    // Se obtienen los indices locales "i" y "j"
    int index_i = obtener_nodos_del_parche[0] % nx;
    int index_j = obtener_nodos_del_parche[0] / nx;

    // Frontera este (\hat{i})
    if (index_i == nx - 1) {
        vecUnitNormal = 1.0;
    }

    // Frontera oeste (-\hat{i})
    if (index_i == 0) {
        vecUnitNormal = -1.0;
    }

    // Frontera norte (\hat{j})
    if (index_j == ny - 1) {
        vecUnitNormal = 1.0;
    }

    // Frontera sur (-\hat{j})
    if (index_j == 0) {
        vecUnitNormal = -1.0;
    }

}

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
)
{
    // Vector 2D que almacena las frontera para evitar repetir codigo
    std::vector<std::vector<Parches_Flujo_de_Masa>> listas;

    listas.emplace_back(parches_norte_FM);
    listas.emplace_back(parches_sur_FM);
    listas.emplace_back(parches_este_FM);
    listas.emplace_back(parches_oeste_FM);

    for (int i = 0; i < static_cast<int>(listas.size()); ++i) {

        // Direccion "x"
        asignar_condiciones_de_frontera_MDot
        (
            listas[i],
            u_star,
            mDotStar_x,
            lista_Dirichlet_x,
            lista_Zero_Neumann_x
        );

        // Direccion "y"
        asignar_condiciones_de_frontera_MDot
        (
            listas[i],
            v_star,
            mDotStar_y,
            lista_Dirichlet_y,
            lista_Zero_Neumann_y
        );

    }

}


void asignar_condiciones_de_frontera_MDot
(
    const std::vector<Parches_Flujo_de_Masa> & parches,
    const std::vector<double>                & vel_star,
    std::vector<double>                      & mDotStar,
    std::vector<CF_MDot<Dirichlet_MDot>>     & lista_Dirichlet,
    std::vector<CF_MDot<Zero_Neumann_MDot>>  & lista_Zero_Neumann
)
{



}
