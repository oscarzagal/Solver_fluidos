//
// Created by oscar on 22/06/25.
//

#include "Campo.hpp"

#include "condiciones_de_frontera.hpp"
#include "ecuaciones_gobernantes.hpp"
#include "malla_por_bloques.hpp"
#include "config_control.hpp"
#include <vector>


namespace Campo {


/*-----------------------------------------------------------------------------
                            Campo Escalar
-----------------------------------------------------------------------------*/

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
     const Malla::Mallador& malla_,
     std::optional<std::reference_wrapper<Ecuaciones_gobernantes::Base>> ecuacion_,
     const std::string& solver_elegido_
    ) :
    phi_new(phi_new_),
    phi_old(phi_old_),
    Parches_norte(Parches_norte_),
    Parches_sur(Parches_sur_),
    Parches_este(Parches_este_),
    Parches_oeste(Parches_oeste_),
    g_dirichlet(g_dirichlet_),
    g_zero_neumann(g_zero_neumann_),
    malla(malla_),
    ecuacion_b(ecuacion_),
    solver_elegido(solver_elegido_)
    {}

    // Constructor que no admite una referencia a "Ecuaciones_gobernantes::Base &"
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
     const Malla::Mallador& malla_,
     const std::string& solver_elegido_
    ) :
    phi_new(phi_new_),
    phi_old(phi_old_),
    Parches_norte(Parches_norte_),
    Parches_sur(Parches_sur_),
    Parches_este(Parches_este_),
    Parches_oeste(Parches_oeste_),
    g_dirichlet(g_dirichlet_),
    g_zero_neumann(g_zero_neumann_),
    malla(malla_),
    solver_elegido(solver_elegido_)
    {}

    std::vector<Condicion_frontera::Dirichlet> Escalar::obtener_parches_dirichlet() {
        return parches_dirichlet;
    }

    std::vector<std::shared_ptr<Condicion_frontera::Base>> Escalar::obtener_parches_dinamicos() {
        return parches_dinamicos;
    }

    void Escalar::construir_condiciones_de_frontera() {

        // Numero de nodos en x
        const int nx = static_cast<int>(malla.obtener_el_numero_de_nodos(Malla::Nodos::nx));

        // Se modifica el estado de "parches_dirichlet" y "parches_dinamicos"
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
        ecuacion_b.ensamblar();

        // Obtencion de los coeficientes agrupados
        A = ecuacion_b.obtener_coeficientes();

        // Numero de nodos en "x" y "y"
        const int nx = static_cast<int>(malla.obtener_el_numero_de_nodos(Malla::Nodos::nx));
        const int ny = static_cast<int>(malla.obtener_el_numero_de_nodos(Malla::Nodos::ny));

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

/*-----------------------------------------------------------------------------
                            Campo Vectorial
-----------------------------------------------------------------------------*/

// Constructor
Vectorial::Vectorial
(
    const Malla::Mallador& malla_,
    const std::vector<double>& Pstar_,
    const std::array<CF_Dirichlet,limite_num_parches> & g_dirichlet_u_,
    const std::array<CF_Zero_Neumann,limite_num_parches> & g_zero_neumann_u_,
    const std::array<CF_Dirichlet,limite_num_parches> & g_dirichlet_v_,
    const std::array<CF_Zero_Neumann,limite_num_parches> & g_zero_neumann_v_,
    const std::string& solver_elegido_
) :
    malla(malla_),
    Pstar(Pstar_),
    g_dirichlet_u(g_dirichlet_u_),
    g_zero_neumann_u(g_zero_neumann_u_),
    g_dirichlet_v(g_dirichlet_v_),
    g_zero_neumann_v(g_zero_neumann_v_),
    solver_elegido(solver_elegido_),
    nx(malla_.obtener_el_numero_de_nodos(Malla::Nodos::nx)),
    ny(malla_.obtener_el_numero_de_nodos(Malla::Nodos::ny)),
    Parches_norte(malla_.obtener_parches(Malla::Frontera::Norte)),
    Parches_sur(malla_.obtener_parches(Malla::Frontera::Sur)),
    Parches_este(malla_.obtener_parches(Malla::Frontera::Este)),
    Parches_oeste(malla_.obtener_parches(Malla::Frontera::Oeste)),
    u_new(nx*ny,0.0),
    v_new(nx*ny,0.0),
    u_old(nx*ny,0.0),
    v_old(nx*ny,0.0),
    ecuacion_momentum
    (
        nu,
        malla_,
        Pstar_,
        u_new,
        v_new
    ),
    // NOTE: constructor sobrecargado
    u
    (
        Parches_norte,
        Parches_sur,
        Parches_este,
        Parches_oeste,
        u_new,
        u_old,
        g_dirichlet_u_,
        g_zero_neumann_u_,
        malla_,
        solver_elegido_
    ),
    v
    (
        Parches_norte,
        Parches_sur,
        Parches_este,
        Parches_oeste,
        v_new,
        v_old,
        g_dirichlet_v_,
        g_zero_neumann_v_,
        malla_,
        solver_elegido_
    )
{}

// TODO: hacer unos getters para devolver a "parches_dirichlet_*" y "parches_dinamicos_*",
// preferiblemente mediante enums

void Vectorial::construir_condiciones_de_frontera() {

    // Condiciones de frontera para la velocidad "u"
    // Se modifica el estado de "parches_dirichlet_u", "parches_dinamicos_u"
    Condicion_frontera::construir_condiciones_de_frontera
    (
        Parches_norte,
        Parches_sur,
        Parches_este,
        Parches_oeste,
        u_new,
        nx,
        u.parches_dirichlet,
        u.parches_dinamicos,
        g_dirichlet_u,   // Array con las condiciones de frontera
        g_zero_neumann_u // Array con las condiciones de frontera
    );

    // Condiciones de frontera para la velocidad "v"
    // Se modifica el estado de "parches_dirichlet_v", "parches_dinamicos_v"
    Condicion_frontera::construir_condiciones_de_frontera
    (
        Parches_norte,
        Parches_sur,
        Parches_este,
        Parches_oeste,
        v_new,
        nx,
        v.parches_dirichlet,
        v.parches_dinamicos,
        g_dirichlet_v,   // Array con las condiciones de frontera
        g_zero_neumann_v // Array con las condiciones de frontera
    );

}

void Vectorial::construir_ecuacion() {

    // Ensamblar la ecuacion gobernante
    ecuacion_momentum.unir_ecuacion();

    // Obtencion de los coeficientes agrupados
    A_u = ecuacion_momentum.obtener_coeficientes_para_u();
    A_v = ecuacion_momentum.obtener_coeficientes_para_v();


    // Asignacion de los coeficientes agrupados A_u, se altera el estado de "u.campo"
    Solver_lineal::asignar
    (
        nx,
        ny,
        u_new,
        u_old,
        A_u,
        solver_elegido,
        u.campo
    );

    // Asignacion de los coeficientes agrupados A_v, se altera el estado de "v.campo"
    Solver_lineal::asignar
    (
        nx,
        ny,
        v_new,
        v_old,
        A_v,
        solver_elegido,
        v.campo
    );

}

void Vectorial::resolver() {
    u.campo->resolver();
    v.campo->resolver();
}

} // Fin namespace Campos
