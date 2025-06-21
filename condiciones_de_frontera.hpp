//
// Created by oscar on 8/06/25.
//

#ifndef CONDICIONES_DE_FRONTERA_HPP
#define CONDICIONES_DE_FRONTERA_HPP

#include <memory>
#include <vector>
#include <string>

#include "malla_por_bloques.hpp"

namespace Condicion_frontera {

    class Base {
    public:

        // Metodo virtual puro
        virtual void aplicar() = 0;

        // Destructor virtual por defecto
        virtual ~Base() = default;

    };

    class Dirichlet {
    public:

        // Constructor
        Dirichlet
        (
            std::vector<double> &,
            std::vector<int> &,
            double
        );

        void aplicar();

    private:

        // Campo
        std::vector<double>& phi;

        std::vector<int>& nodos_del_parche;

        const double valor;
    };

    class Zero_Neumann : public Base {

    public:

        // Constructor
        Zero_Neumann
        (
            std::vector<double>&,
            const std::vector<int>&,
            const std::string&,
            const int&
        );

        void aplicar() override;

    private:

        // Campo
        std::vector<double>& phi;

        const std::vector<int>& nodos_del_parche;

        // Frontera fisica del dominio
        const std::string& frontera_fisica;

        // Numero de nodos en la direccion "x"
        const int& nx;

    };

    class Fabrica_de_CF {
    public:

        // Los argumentos de la funcion deben de ser los mismos para todos los
        // tipos de condiciones de frontera dinamicos
        static std::unique_ptr<Base> crear
        (
            const std::pair<std::string,int> &tipo,
            const Malla::Mallador::Parche& parche,
            std::vector<double>& phi,
            const int& nx
        );

    };

    // Funcion para clasificar los parches en base a su tipo de condicion de
    // frontera
    void asignar_condiciones_de_frontera
    (
        std::vector<Malla::Mallador::Parche>& parches,
        std::vector<double>& phi,
        const int& nx,
        std::vector<Dirichlet>& lista_dirichlet,
        std::vector<std::unique_ptr<Base>>& lista_parches_dinamicos
    );

    // Funcion que retorna el tipo de CF. En caso de no encontrar el parche
    // retorna su nombre
    std::pair<std::string,int> que_tipo_es(const std::string& nombre);


} // Fin namespace Condicion_frontera

#endif //CONDICIONES_DE_FRONTERA_HPP
