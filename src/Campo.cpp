//
// Created by oscar on 22/06/25.
//

#include "Campo.hpp"

#include <iostream>

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
     std::vector<double>& phi_new_,
     std::vector<double>& phi_old_,
     const std::array<CF_Dirichlet,limite_num_parches> & g_dirichlet_,
     const std::array<CF_Zero_Neumann,limite_num_parches> & g_zero_neumann_,
     Malla::Mallador& malla_,
     Ecuaciones_gobernantes::Base& ecuacion_,
     const std::string& solver_elegido_
    ) :
    Parches_norte(Parches_norte_),
    Parches_sur(Parches_sur_),
    Parches_este(Parches_este_),
    Parches_oeste(Parches_oeste_),
    g_dirichlet(g_dirichlet_),
    g_zero_neumann(g_zero_neumann_),
    malla(malla_),
    ecuacion(ecuacion_),
    solver_elegido(solver_elegido_),
    phi_new(phi_new_),
    phi_old(phi_old_)
    {}

    std::vector<Condicion_frontera::Dirichlet> Escalar::obtener_parches_dirichlet() {
        return parches_dirichlet;
    }

    std::vector<std::shared_ptr<Condicion_frontera::Base>> Escalar::obtener_parches_dinamicos() {
        return parches_dinamicos;
    }

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

        // Numero de nodos en "x" y "y"
        const int nx = static_cast<int>(malla.obtener_coordenadas_tmp_x().size());
        const int ny = static_cast<int>(malla.obtener_coordenadas_tmp_x().size());

       // NOTE: al ser nx y ny variables locales y usarlas para construir un pointer
       // a un objeto de tipo SOR (pasando estas dos por REFERENCIA ya que asi lo
       // declare en el constructor de SOR anteriormente) se joden las cosas al salir
       // del alcance de esta funcion (Escalar::construir_ecuacion()), pues nx y ny
       // se destruyeron. Es importante pasarlos como copia cuando se requiera de su
       // uso fuera del alcance de la funcion.


        // Asignacion de los coeficientes agrupados
        Solver_lineal::asignar
        (
            nx,
            ny,
            phi_new,
            phi_old,
            A,
            solver_elegido,
            campo
        );

    }

    void Escalar::resolver() const {

        campo->resolver();

    }


} // Fin namespace Campos
