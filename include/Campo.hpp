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
#include <memory>

typedef std::vector<Malla::Mallador::Parche> almacenar;

namespace Campo {

class Escalar {
public:

    // Constructor para cuando se tengan dos tipos de CF
    Escalar
    (
        const std::vector<Malla::Mallador::Parche> &,
        const std::vector<Malla::Mallador::Parche> &,
        const std::vector<Malla::Mallador::Parche> &,
        const std::vector<Malla::Mallador::Parche> &,
        std::vector<double> &,
        std::vector<double> &,
        const std::array<CF_Dirichlet, limite_num_parches> &,
        const std::array<CF_Zero_Neumann, limite_num_parches> &,
        const Malla::Mallador &,
        Ecuaciones_gobernantes::Base &,
        const std::string&
    );

    // Constructor para admitir una referencia a "Ecuaciones_gobernantes::Momentum &"
    Escalar
    (
        const std::vector<Malla::Mallador::Parche> &,
        const std::vector<Malla::Mallador::Parche> &,
        const std::vector<Malla::Mallador::Parche> &,
        const std::vector<Malla::Mallador::Parche> &,
        std::vector<double> &,
        std::vector<double> &,
        const std::array<CF_Dirichlet, limite_num_parches> &,
        const std::array<CF_Zero_Neumann, limite_num_parches> &,
        const Malla::Mallador &,
        Ecuaciones_gobernantes::Momentum &,
        const std::string&
    );


    void construir_condiciones_de_frontera();

    void construir_ecuacion();

    void resolver() const;

    std::vector<Condicion_frontera::Dirichlet> obtener_parches_dirichlet();

    std::vector<std::shared_ptr<Condicion_frontera::Base>> obtener_parches_dinamicos();

    std::vector<double>& phi_new; // Campo nuevo
    std::vector<double>& phi_old; // Campo viejo

    // Pointer hacia el solver lineal
    std::unique_ptr<Solver_lineal::Base> campo;

private:

    // Variables del constructor
    std::vector<Malla::Mallador::Parche> Parches_norte;
    std::vector<Malla::Mallador::Parche> Parches_sur;
    std::vector<Malla::Mallador::Parche> Parches_este;
    std::vector<Malla::Mallador::Parche> Parches_oeste;
    const std::array<CF_Dirichlet, limite_num_parches> & g_dirichlet;
    const std::array<CF_Zero_Neumann, limite_num_parches> & g_zero_neumann;
    const Malla::Mallador& malla;
    Ecuaciones_gobernantes::Base& ecuacion; // Referencia hacia el tipo de ecuacion
    const std::string& solver_elegido;

    // Variables adicionales
    // Variables que almacenan las condiciones de frontera dependiendo el tipo
    std::vector<Condicion_frontera::Dirichlet> parches_dirichlet;
    std::vector<std::shared_ptr<Condicion_frontera::Base>> parches_dinamicos;


    // Instancia para almacenar los coeficientes agrupados
    Ecuaciones_gobernantes::A_coef A;

};

class Vectorial {
public:

    Vectorial
    (
        const Malla::Mallador &,
        const std::array<CF_Dirichlet, limite_num_parches> &,
        const std::array<CF_Zero_Neumann, limite_num_parches> &,
        const std::array<CF_Dirichlet, limite_num_parches> &,
        const std::array<CF_Zero_Neumann, limite_num_parches> &,
        Ecuaciones_gobernantes::Momentum &, // No hay necesidad de polimorfismo
        const std::string&
    );

    void construir_condiciones_de_frontera();

    void construir_ecuacion();

    // Se aplica el solver lineal a los campos
    void resolver();

private:

    // Variables del constructor
    const Malla::Mallador& malla;
    const std::array<CF_Dirichlet, limite_num_parches> & g_dirichlet_u;
    const std::array<CF_Zero_Neumann, limite_num_parches> & g_zero_neumann_u;
    const std::array<CF_Dirichlet, limite_num_parches> & g_dirichlet_v;
    const std::array<CF_Zero_Neumann, limite_num_parches> & g_zero_neumann_v;
    Ecuaciones_gobernantes::Momentum& ecuacion_momentum;
    const std::string& solver_elegido;

    // Numero de nodos en "x" y "y"
    const int nx, ny;

    // Variables adicionales
    // Parches de frontera
    almacenar Parches_norte;
    almacenar Parches_sur;
    almacenar Parches_este;
    almacenar Parches_oeste;

    // Instancias para almacenar los coeficientes agrupados
    Ecuaciones_gobernantes::A_coef A_u, A_v;


public:

    std::vector<double> u_new, v_new; // Campos nuevos
    std::vector<double> u_old, v_old; // Campos viejos

    Campo::Escalar u; // Campo de velocidad en "x"
    Campo::Escalar v; // Campo de velocidad en "y"

    // Variables que almacenan las condiciones de frontera dependiendo el tipo
    std::vector<Condicion_frontera::Dirichlet> parches_dirichlet_u;
    std::vector<std::shared_ptr<Condicion_frontera::Base>> parches_dinamicos_u;

    std::vector<Condicion_frontera::Dirichlet> parches_dirichlet_v;
    std::vector<std::shared_ptr<Condicion_frontera::Base>> parches_dinamicos_v;




};

} // Fin namespace Campo_Escalar



#endif //CAMPOS_HPP
