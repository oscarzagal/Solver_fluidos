//
// Created by oscar on 22/06/25.
//

#include "ecuaciones_gobernantes.hpp"
#include "esquemas_de_discretizacion.hpp"
#include "config_control.hpp"

namespace Ecuaciones_gobernantes {

    // Lista de inicializacion
    Energia::Energia(Malla::Mallador & malla_) : malla(malla_) {}

    void Energia::ensamblar() {

        const int nx = static_cast<int>(malla.obtener_coordenadas_tmp_x().size());
        const int ny = static_cast<int>(malla.obtener_coordenadas_tmp_y().size());

        A_coef A;

        Esquemas_discretizacion::laplaciano(nx,ny,k,A,malla);

        asignar_matriz(A);

    }

    void Energia::asignar_matriz(const A_coef& A_paso) {
        A = A_paso;
    }


    A_coef Energia::obtener_coeficientes() {
        return A;
    }


} // Fin namespace Ecuaciones_gobernantes