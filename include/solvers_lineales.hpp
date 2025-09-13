//
// Created by oscar on 13/06/25.
//

#ifndef SOLVERS_LINEALES_HPP
#define SOLVERS_LINEALES_HPP

#include <iostream>
#include <memory>
#include <stdexcept>
#include <variant>
#include <vector>
#include "Campo.hpp"

namespace Solver_lineal {

    class SOR {
    public:

        // Constructor
        SOR
        (
            const int nx_,
            const int ny_,
            const double lambda_,
            std::vector<double> &phi_,
            std::vector<double> &phi_old_

        ) :
            nx(nx_),
            ny(ny_),
            lambda(lambda_),
            phi(phi_),
            phi_old(phi_old_)
            {}

            // Se lo pasas por parametro porque ese objeto se obtiene en "runtime"
            void resolver(const Campo::A_coef& A) {

            for (int j=1; j<ny-1; ++j) {
              for (int i=1; i<nx-1; ++i) {
                phi[i+nx*j]=lambda*((phi[i+1+nx*j]*A.ae[i+nx*j]+phi[i-1+nx*j]*A.aw[i+nx*j]
                  +phi[i+nx*(j+1)]*A.an[i+nx*j]+phi[i+nx*(j-1)]*A.as[i+nx*j]+A.b[i+nx*j])
                  /A.ac[i+nx*j])+(1.0-lambda)*phi_old[i+nx*j];
              }
            }

        }

    private:

        // Tamaño de la malla
        const int nx;
        const int ny;

        const double lambda;

        // Campos
        std::vector<double>& phi;
        std::vector<double>& phi_old;

    };


    /*-----------------------------------------------------------------------------
            Eleccion del solver lineal (auqnue por el momento solo hay uno)
    -----------------------------------------------------------------------------*/

    // Declaracion del variant
    using solverVariant = std::variant<SOR>;

    inline solverVariant solverElegido
    (
        const int nx,
        const int ny,
        const double lambda,
        std::vector<double>& phi,
        std::vector<double>& phi_old,
        const std::string& solver_elegido
    )
    {
        if (solver_elegido == "SOR") {
            return SOR(nx, ny, lambda, phi, phi_old);
        } else {
            throw std::runtime_error("Solver no soportado: " + solver_elegido);
        }

    }

    /*-----------------------------------------------------------------------------
          Fin Eleccion del solver lineal (auqnue por el momento solo hay uno)
    -----------------------------------------------------------------------------*/


} // Fin namespace Solver_lineal

#endif //SOLVERS_LINEALES_HPP
