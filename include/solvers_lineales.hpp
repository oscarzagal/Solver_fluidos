//
// Created by oscar on 13/06/25.
//

#ifndef SOLVERS_LINEALES_HPP
#define SOLVERS_LINEALES_HPP

#include <iostream>
#include <memory>
#include <vector>
#include <stdexcept>

#include "config_control.hpp"
#include "ecuaciones_gobernantes.hpp"

namespace Solver_lineal {

    class Base {
    public:

        // Metodo virtual puro
        virtual void resolver() = 0;

        // Destructor virtual por defecto
        virtual ~Base() = default;

    };

    class SOR : public Base {
    public:

        // Constructor
        SOR
        (
            const int nx_,
            const int ny_,
            std::vector<double> &phi_,
            std::vector<double> &phi_old_,
            const Ecuaciones_gobernantes::A_coef& A_
        ) :
            nx(nx_),
            ny(ny_),
            phi(phi_),
            phi_old(phi_old_),
            A(A_) {

        }

        void resolver() override {

            for (int j=1; j<ny-1; ++j) {
              for (int i=1; i<nx-1; ++i) {
                phi[i+nx*j]=lambda_T*((phi[i+1+nx*j]*A.ae[i+nx*j]+phi[i-1+nx*j]*A.aw[i+nx*j]
                  +phi[i+nx*(j+1)]*A.an[i+nx*j]+phi[i+nx*(j-1)]*A.as[i+nx*j]+A.b[i+nx*j])
                  /A.ap[i+nx*j])+(1.0-lambda_T)*phi_old[i+nx*j];
              }
            }

        }

    private:

        // Tamaño de la malla
        const int nx;
        const int ny;

        // Campos
        std::vector<double>& phi;
        std::vector<double>& phi_old;

        // Coeficientes agrupados
        const Ecuaciones_gobernantes::A_coef& A;

    };

    inline void asignar
    (
        const int &nx,
        const int &ny,
        std::vector<double>& phi,
        std::vector<double> &phi_old,
        const Ecuaciones_gobernantes::A_coef& A,
        const std::string& solver_elegido,
        std::unique_ptr<Base>& campo
    )
    {

        if (solver_elegido == "SOR") {
            campo = std::make_unique<SOR>(nx,ny,phi,phi_old,A);
        } else {
            throw std::runtime_error("No existe el solver " + solver_elegido);
        }

    }


} // Fin namespace Solver_lineal

#endif //SOLVERS_LINEALES_HPP
