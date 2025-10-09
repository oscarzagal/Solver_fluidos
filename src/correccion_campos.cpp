//
// Created by oscar on 04/10/25.
//

#include "correccion_campos.hpp"
#include "config_control.hpp"
#include "utilidades.hpp"
#include "variables_discretizacion.hpp"
#include "esquemas_de_discretizacion.hpp"
#include <algorithm>
#include <execution>

// Constructor
Correccion::Correccion
(
    const Malla::Mallador & malla_,
    const Coeficiente_d   & coef_d_,
    Campo::MDotStar       & mdotstar_,
    Campo::Momentum       & velU_,
    Campo::Presion        & presion_
 ) :
    malla(malla_),
    coef_d(coef_d_),
    mdotstar(mdotstar_),
    velU(velU_),
    presion(presion_),
    nx(malla_.obtener_el_numero_de_nodos(Malla::Nodos::nx)),
    ny(malla_.obtener_el_numero_de_nodos(Malla::Nodos::ny)),
    inter(Malla::Mallador::obtener_factores_de_interpolacion(malla_)),
    gradPrime_times_vol(nx, ny),
    vol(malla_.obtener_volumenes()),
    celdas_interiores((nx - 2) * (ny - 2), 0)

{}


/*-----------------------------------------------------------------------------
                            Funciones miembro
-----------------------------------------------------------------------------*/

void Correccion::obtener_celdas_interiores() {

    int iter = 0;
    for (int j = 1 ; j < ny - 1 ; ++j) {
      for (int i = 1 ; i < nx - 1 ; ++i) {
          celdas_interiores[iter] = idx(i, j, nx);
          iter++;
      }
    }

}


void Correccion::corregir() {

    Discretizacion::Explicita::gradiente(nx, ny, inter, gradPrime_times_vol, presion.Pprime, malla);

    // Coordenadas de los centroides de las celdas
    const auto& x = malla.obtener_coord_pers_x();
    const auto& y = malla.obtener_coord_pers_y();

    // Recordar que las deltas representan el area de las caras de los elementos
    // computacionales
    const auto& deltay = malla.deltay;
    const auto& deltax = malla.deltax;

    auto& u_star = velU.u_star;
    auto& v_star = velU.v_star;

    const auto& dC_u = coef_d.dC_u;
    const auto& dE_u = coef_d.dE_u;
    const auto& dC_v = coef_d.dC_v;
    const auto& dN_v = coef_d.dN_v;

    /*-----------------------------------------------------------------------------
                           Inicio Lambda iterativa
    -----------------------------------------------------------------------------*/

    auto bucle = [&](int Centro) {

        const int Este  = Centro + 1;
        const int Norte = Centro + nx;

        const int i = Centro % nx;
        const int j = Centro / nx;

        // Area de la caras locales "e" y "n"
        const double S_e = deltay[j];
        const double S_n = deltax[i];

        // Gradientes de presion de correccion en la celda "C"
        const double gradPrime_Cx = gradPrime_times_vol.grad_x_vol[Centro] / vol[Centro];
        const double gradPrime_Cy = gradPrime_times_vol.grad_y_vol[Centro] / vol[Centro];

        // Factores de interpolacion
        const double gx = inter.ge[Centro];
        const double gy = inter.gn[Centro];

        // Coeficiente "d" interpolado
        const double d_interp_x = interpolar(dC_u[Centro], dE_u[Centro], gx);
        const double d_interp_y = interpolar(dC_v[Centro], dN_v[Centro], gy);

        // Aliases para el flujo de masa
        double& mdotstar_e = mdotstar.mDotStar_x[Este];
        double& mdotstar_n = mdotstar.mDotStar_y[Norte];

        // Presiones de correccion
        const double Pprime_C = presion.Pprime[Centro];
        const double Pprime_E = presion.Pprime[Este];
        const double Pprime_N = presion.Pprime[Norte];

        // δx_{CE} y δX_{CN}
        const double delta_x_CE = x[Este] - x[Centro];
        const double delta_y_CN = y[Norte] - y[Centro];

        // Gradiente de presion sobre las caras "e" y "n"
        const double gradPprime_e = (Pprime_E - Pprime_C) / delta_x_CE;
        const double gradPprime_n = (Pprime_N - Pprime_C) / delta_y_CN;

        // Alias para la presion
        double& Pstar_C = presion.P_star[Centro];

        /*-----------------------------------------------------------------------------
                                    Correccion de campos
        -----------------------------------------------------------------------------*/

        // Velocidad en los centroides de los elementos
        u_star[Centro] = u_star[Centro] - dC_u[Centro] * gradPrime_Cx;
        v_star[Centro] = v_star[Centro] - dC_v[Centro] * gradPrime_Cy;

        // Flujos de masa
        mdotstar_e = mdotstar_e - d_interp_x * gradPprime_e * S_e;
        mdotstar_n = mdotstar_n - d_interp_y * gradPprime_n * S_n;

        // Presion
        Pstar_C = Pstar_C + lambda_P * Pprime_C;

        /*-----------------------------------------------------------------------------
                                    Fin Correccion de campos
        -----------------------------------------------------------------------------*/

    };

    /*-----------------------------------------------------------------------------
                                Fin Lambda iterativa
    -----------------------------------------------------------------------------*/

    // DEBUG
    // std::cout << "celdas_interiores.begin() = " << *celdas_interiores.begin() << "\n";
    // std::cout << "celdas_interiores.end() - 1 = " << *(celdas_interiores.end() - 1) << "\n";

    // Paralelizacion con "for_each"
    std::for_each
    (
#ifdef PARALLEL
        std::execution::par,
#else
        std::execution::seq,
#endif
        celdas_interiores.begin(),
        celdas_interiores.end(),
        bucle
    );



    /*-----------------------------------------------------------------------------
                        Actualizar condiciones de frontera
    -----------------------------------------------------------------------------*/

    // Actualizacion de los parches dinamicos para la presion dinamica
    for (int i = 0; i < static_cast<int>(presion.lista_parches_dinamicos.size()); ++i) {
        presion.lista_parches_dinamicos[i]->aplicar();
    }

    /*-----------------------------------------------------------------------------
                      Fin Actualizar condiciones de frontera
    -----------------------------------------------------------------------------*/

}


/*-----------------------------------------------------------------------------
                            Fin Funciones miembro
-----------------------------------------------------------------------------*/
