//
// Created by oscar on 19/06/25.
//

#include "solvers_lineales.hpp"

#include <stdexcept>

#include "config_control.hpp"

namespace Solver_lineal {

    SOR::SOR
    (
        const int & nx_,
        const int & ny_,
        std::vector<double> & phi_,
        const std::vector<double> & phi_old_,
        const std::vector<double> & ap_,
        const std::vector<double> & ae_,
        const std::vector<double> & aw_,
        const std::vector<double> & an_,
        const std::vector<double> & as_,
        const std::vector<double> & b_
    ) :
        nx(nx_),
        ny(ny_),
        phi(phi_),
        phi_old(phi_old_),
        ap(ap_),
        ae(ae_),
        aw(aw_),
        an(an_),
        as(as_),
        b(b_)
    {}

    void SOR::resolver() {

       for (int j=1; j<ny-1; ++j) {
         for (int i=1; i<nx-1; ++i) {
           phi[i+nx*j]=lambdaT*((phi[i+1+nx*j]*ae[i+nx*j]+phi[i-1+nx*j]*aw[i+nx*j]
             +phi[i+nx*(j+1)]*an[i+nx*j]+phi[i+nx*(j-1)]*as[i+nx*j]+b[i+nx*j])
             /ap[i+nx*j])+(1.0-lambdaT)*phi_old[i+nx*j];
         }
       }

    }

    std::unique_ptr<Base> Fabrica_de_solvers::crear
    (
        const int &nx,
        const int &ny,
        std::vector<double> &phi,
        const std::vector<double> &phi_old,
        const std::vector<double> &ap,
        const std::vector<double> &ae,
        const std::vector<double> &aw,
        const std::vector<double> &an,
        const std::vector<double> &as,
        const std::vector<double> &b
    )
    {
        if (solver_elegido == "SOR") {
            return std::make_unique<SOR>(nx,ny,phi,phi_old,ap,ae,aw,an,as,b);
        }

        return {};
    }


    void asignar
    (
        const int &nx,
        const int &ny,
        std::vector<double>& phi,
        const std::vector<double> &phi_old,
        const std::vector<double> &ap,
        const std::vector<double> &ae,
        const std::vector<double> &aw,
        const std::vector<double> &an,
        const std::vector<double> &as,
        const std::vector<double> &b,
        std::unique_ptr<Base>& campo
    )
    {
        if (solver_elegido == "SOR") {
            campo = Fabrica_de_solvers::crear(nx,ny,phi,phi_old,ap,ae,aw,an,as,b);
        } else {
            throw std::runtime_error("No existe el solver " + solver_elegido);
        }

    }


} // Fin namespace Solver_lineal