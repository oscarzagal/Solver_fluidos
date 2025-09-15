//
// Created by oscar on 1/06/25.
//


#include "malla_por_bloques.hpp"
#include "config_malla.hpp"
#include "condiciones_de_frontera.hpp"
#include "condiciones_de_frontera_MDot.hpp"
#include "Campo.hpp"
#include "variables_discretizacion.hpp"
#include "ecuacion_momentum.hpp"
#include "config_control.hpp"
#include "solvers_lineales.hpp"
#include "calculo_del_error.hpp"
#include "escritura.hpp"
#include "config_CF.hpp"
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <vector>

constexpr bool debug  = false; // Parches flujo de masa
constexpr bool debug2 = false; // Listas para las especializaciones y flujo de masa

int main() {

    //--------------------------Creacion de la malla----------------------------

    // Creacion del objeto que representa la malla actual
    Malla::Mallador malla
    (
        nodos_en_x,
        nodos_en_y,
        coordenadas_en_x,
        coordenadas_en_y,
        nombres_frontera_norte,
        nombres_frontera_sur,
        nombres_frontera_este,
        nombres_frontera_oeste
    );

    // Primero tienes que preparar para luego obtener las coordenadas persistentes
    malla.preparar_coordenadas_persistentes();

    // Numero de nodos en x
    const int nx = malla.obtener_el_numero_de_nodos(Malla::Nodos::nx);

    // Numero de nodos en y
    const int ny = malla.obtener_el_numero_de_nodos(Malla::Nodos::ny);

    // for (int i = 0; i < ny; ++i) {
    //     printf("deltay[%d] = %f\n", i, malla.deltay[i]);
    // }

    // std::cout << "ge en el main.cpp \n";
    // for (int j = 0; j < ny; ++j) {
    //     for (int i = 0; i < nx; ++i) {
    //         printf("ge[%d] = %f\n", i + nx * j, inter.ge[i + nx * j]);
    //     }
    // }

    // Coordenadas persistentes
    const std::vector<double> x = malla.obtener_coord_pers_x();
    const std::vector<double> y = malla.obtener_coord_pers_y();

    // Solo se debe usar una vez para evitar duplicaciones y cosas raras
    malla.preparar_parches_fronteras();

    // Obtener parches para su uso futuro
    almacenar Parches_norte = malla.obtener_parches(Malla::Frontera::Norte);
    almacenar Parches_sur   = malla.obtener_parches(Malla::Frontera::Sur);
    almacenar Parches_este  = malla.obtener_parches(Malla::Frontera::Este);
    almacenar Parches_oeste = malla.obtener_parches(Malla::Frontera::Oeste);

    //------------------------Fin creacion de la malla--------------------------

    /*-----------------------------------------------------------------------------
                            Parches para el flujo de masa
    -----------------------------------------------------------------------------*/

    // TODO: hacer toda esta madre una funcion

    std::vector<Parches_Flujo_de_Masa> parches_norte_FM;
    std::vector<Parches_Flujo_de_Masa> parches_sur_FM;
    std::vector<Parches_Flujo_de_Masa> parches_este_FM;
    std::vector<Parches_Flujo_de_Masa> parches_oeste_FM;

                            // NOTE: Espacio, lo que lleva cada espacio (objeto)
    parches_norte_FM.resize(Parches_norte.size(), Parches_Flujo_de_Masa(nx, ny));
    parches_sur_FM.resize(Parches_sur.size(), Parches_Flujo_de_Masa(nx, ny));
    parches_este_FM.resize(Parches_este.size(), Parches_Flujo_de_Masa(nx, ny));
    parches_oeste_FM.resize(Parches_oeste.size(), Parches_Flujo_de_Masa(nx, ny));

    // NOTE: se debe de cuidar mucho el orden, lo cual hace que el codigo sea muy
    // fragil. Es necesario corregir este defecto en el futuro, por ahora funciona
    // y no es la prioridad.
    for (int i = 0; i < static_cast<int>(parches_norte_FM.size()); ++i) {
        parches_norte_FM[i].obtener_nodos_del_parche = Parches_norte[i].obtener_nodos_del_parche;
        parches_norte_FM[i].obtener_nombre = Parches_norte[i].obtener_nombre;
        parches_norte_FM[i].cortar_nodos_esquina();
        parches_norte_FM[i].calcular_vector_normal_unitario();
        parches_norte_FM[i].asignar_deltas(malla);
        parches_norte_FM[i].calcular_desfase();
        parches_norte_FM[i].tipo_de_CF = parches_norte_FM[i].añadir_tipo_de_CF(g_dirichlet_u, g_zero_neumann_u);
    }

    for (int i = 0; i < static_cast<int>(parches_sur_FM.size()); ++i) {
        parches_sur_FM[i].obtener_nodos_del_parche = Parches_sur[i].obtener_nodos_del_parche;
        parches_sur_FM[i].obtener_nombre = Parches_sur[i].obtener_nombre;
        parches_sur_FM[i].cortar_nodos_esquina();
        parches_sur_FM[i].calcular_vector_normal_unitario();
        parches_sur_FM[i].asignar_deltas(malla);
        parches_sur_FM[i].calcular_desfase();
        parches_sur_FM[i].tipo_de_CF = parches_sur_FM[i].añadir_tipo_de_CF(g_dirichlet_u, g_zero_neumann_u);
    }

    for (int i = 0; i < static_cast<int>(parches_este_FM.size()); ++i) {
        parches_este_FM[i].obtener_nodos_del_parche = Parches_este[i].obtener_nodos_del_parche;
        parches_este_FM[i].obtener_nombre = Parches_este[i].obtener_nombre;
        parches_este_FM[i].cortar_nodos_esquina();
        parches_este_FM[i].calcular_vector_normal_unitario();
        parches_este_FM[i].asignar_deltas(malla);
        parches_este_FM[i].calcular_desfase();
        parches_este_FM[i].tipo_de_CF = parches_este_FM[i].añadir_tipo_de_CF(g_dirichlet_u, g_zero_neumann_u);
    }

    for (int i = 0; i < static_cast<int>(parches_oeste_FM.size()); ++i) {
        parches_oeste_FM[i].obtener_nodos_del_parche = Parches_oeste[i].obtener_nodos_del_parche;
        parches_oeste_FM[i].obtener_nombre = Parches_oeste[i].obtener_nombre;
        parches_oeste_FM[i].cortar_nodos_esquina();
        parches_oeste_FM[i].calcular_vector_normal_unitario();
        parches_oeste_FM[i].asignar_deltas(malla);
        parches_oeste_FM[i].calcular_desfase();
        parches_oeste_FM[i].tipo_de_CF = parches_oeste_FM[i].añadir_tipo_de_CF(g_dirichlet_u, g_zero_neumann_u);
    }


    if (debug) {
        constexpr int numero_parche = 1;

#define FRONTERA_PARCHE parches_este_FM

        std::cout << "Frontera fisica: " << FRONTERA_PARCHE[numero_parche].frontera_fisica << "\n";
        std::cout << "Vector unitario: " << FRONTERA_PARCHE[numero_parche].vecUnitNormal << "\n";
        std::cout << "Desfase del parche: " << FRONTERA_PARCHE[numero_parche].desfase << "\n";
        std::cout << "Tipo de CF: " << FRONTERA_PARCHE[numero_parche].tipo_de_CF << "\n";
        std::cout << "Nodos del parche " << FRONTERA_PARCHE[numero_parche].obtener_nombre << " (copia) luego de cortar: \n";
        for (const int nodo : FRONTERA_PARCHE[numero_parche].obtener_nodos_del_parche) {
            std::cout << nodo << " ";
        }

        std::cout << "\n\n";

        std::cout << "Delta asignada: " << FRONTERA_PARCHE[numero_parche].que_delta_fue_asignada << "\n";
        std::cout << "Deltas: \n";
        for (int i = 0; i < (int)FRONTERA_PARCHE[numero_parche].delta.size(); ++i) {
            std::cout << FRONTERA_PARCHE[numero_parche].delta[i] << " ";
        }

    } // Fin debug



    /*-----------------------------------------------------------------------------
                            Parches para el flujo de masa
    -----------------------------------------------------------------------------*/



    /*-----------------------------------------------------------------------------
                            Inicializacion de campos
    -----------------------------------------------------------------------------*/

    Campo::Momentum velU(nx, ny, 5.2, 1.0);
    Campo::Presion  presion(nx, ny, 0.0);
    Campo::velFace  velface(nx, ny, 0.0);
    Campo::MDotStar mdotstar(nx, ny, 0.0);

    Gradiente grad(nx, ny);               // Gradiente de presion por volumen
    fluxes_difusivos flux_dif(nx, ny);
    fluxes_convectivos flux_conv(nx, ny); // Implica el calculo del flujo de masa


    /*-----------------------------------------------------------------------------
                           Fin inicializacion de campos
    -----------------------------------------------------------------------------*/



    /*-----------------------------------------------------------------------------
                      Asignacion de condiciones de frontera
    -----------------------------------------------------------------------------*/

    // Construccion de las condiciones de frontera para la velocidad "u" y "v"
    // NOTE: se podria usar una funcion contenedor
    construir_condiciones_de_frontera
    (
        Parches_norte,
        Parches_sur,
        Parches_este,
        Parches_oeste,
        velU.u_star,
        nx,
        velU.lista_parches_dirichlet_u,
        velU.lista_parches_dinamicos_u,
        g_dirichlet_u,
        g_zero_neumann_u
    );

    construir_condiciones_de_frontera
    (
        Parches_norte,
        Parches_sur,
        Parches_este,
        Parches_oeste,
        velU.v_star,
        nx,
        velU.lista_parches_dirichlet_v,
        velU.lista_parches_dinamicos_v,
        g_dirichlet_v,
        g_zero_neumann_v
    );


    // Construccion de las condiciones de frontera para la presion "P"
    construir_condiciones_de_frontera
    (
        Parches_norte,
        Parches_sur,
        Parches_este,
        Parches_oeste,
        presion.P_star,
        nx,
        presion.lista_parches_dirichlet,
        presion.lista_parches_dinamicos,
        g_dirichlet_P,
        g_zero_neumann_P
    );

    // Construccion de las condiciones de frontera para el flujo de masa.
    construir_CF_flujo_de_masa
    (
        parches_norte_FM,
        parches_sur_FM,
        parches_este_FM,
        parches_oeste_FM,
        velU.u_star,
        velU.v_star,
        mdotstar.mDotStar_x,
        mdotstar.mDotStar_y,
        mdotstar.lista_Dirichlet_x,
        mdotstar.lista_Zero_Neumann_x,
        mdotstar.lista_Dirichlet_y,
        mdotstar.lista_Zero_Neumann_y
    );


    if (debug2)
    {
        std::cout << "\n\n";
        std::cout << "Tamaño de lista_Dirichlet_x: " << mdotstar.lista_Dirichlet_x.size() << "\n";
        std::cout << "Tamaño de lista_Zero_Neumann_x: " << mdotstar.lista_Zero_Neumann_x.size() << "\n";

        std::cout << "Tamaño de lista_Dirichlet_y: " << mdotstar.lista_Dirichlet_y.size() << "\n";
        std::cout << "Tamaño de lista_Zero_Neumann_y: " << mdotstar.lista_Zero_Neumann_y.size() << "\n";


        for (int j = 0 ; j < ny ; ++j) {
            for (int i = 0 ; i < nx ; ++i) {
                printf("mDotStar_x[%d] = %f \n", i + nx * j, mdotstar.mDotStar_x[i + nx * j]);
                // printf("mDotStar_y[%d] = %f \n", i + nx * j, mdotstar.mDotStar_y[i + nx * j]);
            }
        }

    }

    mdotstar.lista_Dirichlet_x[0].aplicar();
    mdotstar.lista_Dirichlet_x[1].aplicar();
    mdotstar.lista_Zero_Neumann_x[0].aplicar();
    mdotstar.lista_Zero_Neumann_x[1].aplicar();

    mdotstar.lista_Dirichlet_y[0].aplicar();
    mdotstar.lista_Dirichlet_y[1].aplicar();
    mdotstar.lista_Zero_Neumann_y[0].aplicar();
    mdotstar.lista_Zero_Neumann_y[1].aplicar();


    /*-----------------------------------------------------------------------------
                      Fin asignacion condiciones de frontera
    -----------------------------------------------------------------------------*/



    /*-----------------------------------------------------------------------------
                               Armado de ecuaciones
    -----------------------------------------------------------------------------*/

    // Instancia de la ecuacion de momentum
    Ecuacion_Momentum ecuacion_momentum(malla, velU, presion, mdotstar, grad, flux_dif, flux_conv);

    ecuacion_momentum.calcular_conductancia_difusiva();

    ecuacion_momentum.resolver();

    for (int j = 0 ; j < ny ; ++j) {
      for (int i = 0 ; i < nx ; ++i) {
          printf("vel_u[%d] = %f\n", i + nx * j, velU.u_star[i + nx * j]);
      }
    }

    /*-----------------------------------------------------------------------------
                            Fin Armado de ecuaciones
    -----------------------------------------------------------------------------*/



    /*-----------------------------------------------------------------------------
                               Inicio bucle SIMPLE
    -----------------------------------------------------------------------------*/



    /*-----------------------------------------------------------------------------
                               Fin Inicio bucle SIMPLE
    -----------------------------------------------------------------------------*/




    return 0;
}
