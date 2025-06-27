//
// Created by oscar on 22/06/25.
//

#include "Campo.hpp"
#include "condiciones_de_frontera.hpp"
#include "malla_por_bloques.hpp"

namespace Campo {

    // Constructor
    Escalar::Escalar
    (
     const std::vector<Malla::Mallador::Parche> & Parches_norte_,
     const std::vector<Malla::Mallador::Parche> & Parches_sur_,
     const std::vector<Malla::Mallador::Parche> & Parches_este_,
     const std::vector<Malla::Mallador::Parche> & Parches_oeste_,
     const std::array<CF_Dirichlet,limite_num_parches> & g_dirichlet_,
     const std::array<CF_Zero_Neumann,limite_num_parches> & g_zero_neumann_,
     Malla::Mallador& malla_,
     Ecuaciones_gobernantes::Base& ecuacion_
    ) :
    Parches_norte(Parches_norte_),
    Parches_sur(Parches_sur_),
    Parches_este(Parches_este_),
    Parches_oeste(Parches_oeste_),
    g_dirichlet(g_dirichlet_),
    g_zero_neumann(g_zero_neumann_),
    malla(malla_),
    ecuacion(ecuacion_)
    {}

    void Escalar::construir_condiciones_de_frontera() {

        // Numero de nodos en x
        const int nx = static_cast<int>(malla.obtener_coordenadas_tmp_x().size());

        Condicion_frontera::construir_condiciones_de_frontera
        (
            Parches_norte,
            Parches_sur,
            Parches_este,
            Parches_oeste,
            phi_new,
            nx,
            parches_dirichlet,
            parches_dinamicos,
            g_dirichlet,
            g_zero_neumann
        );

    }

    void Escalar::construir_ecuacion() {

        // Ensamblar la ecuacion gobernante
        ecuacion.ensamblar();

        // Obtencion de los coeficientes agrupados
        A = ecuacion.obtener_coeficientes();

    }

    void Escalar::resolver() {

    }






} // Fin namespace Campos
