//
// Created by oscar on 1/06/25.
//

#include <iostream>
#include <string>

#include "malla_por_bloques.hpp"
#include "config_malla.hpp"
#include "condiciones_de_frontera.hpp"
#include "config_control.hpp"
#include "solvers_lineales.hpp"
#include "calculo_del_error.hpp"
#include "escritura.hpp"
#include "config_CF.hpp"
#include "Campo.hpp"
#include "ecuaciones_gobernantes.hpp"


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
    const int nx = static_cast<int>(malla.obtener_coordenadas_tmp_x().size());

    // Numero de nodos en y
    const int ny = static_cast<int>(malla.obtener_coordenadas_tmp_y().size());

    // Coordenadas persistentes
    const std::vector<double> x = malla.obtener_coord_pers_x();
    const std::vector<double> y = malla.obtener_coord_pers_y();

    // Solo se debe usar una vez para evitar duplicaciones y cosas raras
    malla.preparar_parches_fronteras();

    // Obtener parches para su uso futuro
    almacenar Parches_norte=malla.obtener_parches(Malla::Frontera::Norte);
    almacenar Parches_sur=malla.obtener_parches(Malla::Frontera::Sur);
    almacenar Parches_este=malla.obtener_parches(Malla::Frontera::Este);
    almacenar Parches_oeste=malla.obtener_parches(Malla::Frontera::Oeste);

    //------------------------Fin creacion de la malla--------------------------

    //--------------------------Definicion de campos----------------------------

    // Campos de temperaturas
    std::vector<double> T(nx*ny,0.0);
    std::vector<double> Told(nx*ny,0.0);

    // Otra inicializacion de prueba
    for (int j = 1; j < ny-1; ++j) {
      for (int i = 1; i < nx-1; ++i) {
          T[i+nx*j]=24.0;
      }
    }

    //------------------------Fin definicion de campos--------------------------

    //-----------------Asignacion de condiciones de frontera--------------------

    // Instancia de la ecuacion de energia que sirve para ensamblar la ecuacion
    Ecuaciones_gobernantes::Energia ecuacion_energia(malla);

    // Creacion del campo escalar T
    Campo::Escalar Temp
    (
        Parches_norte,
        Parches_sur,
        Parches_este,
        Parches_oeste,
        T,
        Told,
        g_dirichlet_T,
        g_zero_neumann_T,
        malla,
        ecuacion_energia,
        solver_elegido_Temp
    );


    Temp.construir_condiciones_de_frontera();

    std::vector<Condicion_frontera::Dirichlet> lista_parches_dirichlet_T
    = Temp.obtener_parches_dirichlet();

    std::vector<std::shared_ptr<Condicion_frontera::Base>> lista_parches_dinamicos_T
    = Temp.obtener_parches_dinamicos();

    Temp.construir_ecuacion();

    //--------------Fin de asignacion de condiciones de frontera----------------

    // Aplicacion de las condiciones de frontera al campo

    for (int i=0; i<lista_parches_dirichlet_T.size(); ++i ) {
        lista_parches_dirichlet_T[i].aplicar();
    }

    int numit=0;
    double error_mayor=1.0;

    while (error_mayor>tolerancia) {

        Temp.resolver();

        // Actualizacion de condiciones de frontera dinamicas
        for (int i = 0; i < lista_parches_dinamicos_T.size(); ++i) {
            lista_parches_dinamicos_T[i]->aplicar();
        }

        error_mayor=calcular_error_mayor(nx,ny,T,Told);

        ++numit;

        printf("Iteracion = %d, Error maximo = %e\n",numit,error_mayor);

        Told=T;

        if (numit==num_iteraciones_max) break;

    }

    std::cout << "\n\n";

    std::cout << "Numero de elementos en x: " << nx << "\n";
    std::cout << "Numero de elementos en y: " << ny << "\n";

    escribir("T.dat","T",x,y,nx,ny,T);

    return 0;
}
