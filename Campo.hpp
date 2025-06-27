//
// Created by oscar on 22/06/25.
//

#ifndef CAMPOS_HPP
#define CAMPOS_HPP


#include "condiciones_de_frontera.hpp"
#include "malla_por_bloques.hpp"
#include "config_CF.hpp"
#include "solvers_lineales.hpp"
#include "ecuaciones_gobernantes.hpp"

namespace Campo {

    class Escalar {
    public:
        // Constructor
        Escalar
        (
            const std::vector<Malla::Mallador::Parche> &,
            const std::vector<Malla::Mallador::Parche> &,
            const std::vector<Malla::Mallador::Parche> &,
            const std::vector<Malla::Mallador::Parche> &,
            const std::array<CF_Dirichlet, limite_num_parches> &,
            const std::array<CF_Zero_Neumann, limite_num_parches> &,
            Malla::Mallador &,
            Ecuaciones_gobernantes::Base &
        );

        void construir_condiciones_de_frontera();

        void construir_ecuacion();

        void resolver();

        std::vector<double> obtener_campo();

        std::vector<Condicion_frontera::Dirichlet> obtener_parches_dirichlet();

        std::vector<Condicion_frontera::Base> obtener_parches_dinamicos();

    private:

        // Variables del constructor
        std::vector<Malla::Mallador::Parche> Parches_norte;
        std::vector<Malla::Mallador::Parche> Parches_sur;
        std::vector<Malla::Mallador::Parche> Parches_este;
        std::vector<Malla::Mallador::Parche> Parches_oeste;
        const std::array<CF_Dirichlet, limite_num_parches> & g_dirichlet;
        const std::array<CF_Zero_Neumann, limite_num_parches> & g_zero_neumann;
        Malla::Mallador& malla;
        Ecuaciones_gobernantes::Base& ecuacion; // Pointer hacia el tipo de ecuacion


        // Variables adicionales
        std::vector<double> phi_new,phi_old; // Campo nuevo y campo viejo
        std::vector<Condicion_frontera::Dirichlet> parches_dirichlet;
        std::vector<std::unique_ptr<Condicion_frontera::Base>> parches_dinamicos;

        // Pointer hacia el solver lineal
        std::unique_ptr<Solver_lineal::Base> campo;

        // Instancia a "Energia::A_energia" para almacenar los coeficientes agrupados
        Ecuaciones_gobernantes::A_coef A;


    };

} // Fin namespace Campo_Escalar



#endif //CAMPOS_HPP
