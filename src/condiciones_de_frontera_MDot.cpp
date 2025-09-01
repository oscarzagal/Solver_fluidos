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
