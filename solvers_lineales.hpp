//
// Created by oscar on 13/06/25.
//

#ifndef SOLVERS_LINEALES_HPP
#define SOLVERS_LINEALES_HPP

#include <memory>
#include <vector>

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
            const int&,
            const int&,
            std::vector<double>&,
            const std::vector<double>&,
            const std::vector<double>&,
            const std::vector<double>&,
            const std::vector<double>&,
            const std::vector<double>&,
            const std::vector<double>&,
            const std::vector<double>&
        );

        void resolver() override;

    private:

        // Tamaño de la malla
        const int& nx;
        const int& ny;

        // Campos
        std::vector<double>& phi;
        const std::vector<double>& phi_old;

        // Coeficientes agrupados
        const std::vector<double>& ap;
        const std::vector<double>& ae;
        const std::vector<double>& aw;
        const std::vector<double>& an;
        const std::vector<double>& as;
        const std::vector<double>& b;

    };

    class Fabrica_de_solvers {
    public:
        static std::unique_ptr<Base> crear
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
            const std::vector<double> &b
        );
    };


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
    );


} // Fin namespace Solver_lineal

#endif //SOLVERS_LINEALES_HPP
