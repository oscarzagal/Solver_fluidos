//
// Created by oscar on 8/06/25.
//

#ifndef CONDICIONES_DE_FRONTERA_HPP
#define CONDICIONES_DE_FRONTERA_HPP

#include <memory>
#include <vector>
#include <string>

#include "malla_por_bloques.hpp"
#include "config_CF.hpp"

namespace Condicion_frontera {

// TODO: Implementar TEMPLATE SPECIALIZATION para las condiciones de frontera
// luego de ver que el codigo funciona. Como esto funciona NO SE TOCA por el
// momento.

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
            std::vector<int> &,
            double,
            std::vector<double> &
        );

        void aplicar();

        std::vector<int>& nodos_del_parche;

        const double valor;

    private:

        // Campo
        std::vector<double>& phi;

    };

    class Zero_Neumann : public Base {

    public:

        // Constructor
        Zero_Neumann
        (
            std::vector<double>&,
            const std::vector<int>&,
            const std::string&,
            int
        );

        void aplicar() override;

    private:

        // Campo
        std::vector<double>& phi;

        const std::vector<int>& nodos_del_parche;

        // Frontera fisica del dominio
        const std::string& frontera_fisica;

        // Numero de nodos en la direccion "x"
        const int nx;

    };

    class Fabrica_de_CF {
    public:

        // Los argumentos de la funcion deben de ser los mismos para todos los
        // tipos de condiciones de frontera dinamicos
        static std::shared_ptr<Base> crear
        (
            const std::pair<std::string,int> &tipo,
            const Malla::Mallador::Parche& parche,
            std::vector<double>& phi,
            int nx,
            const std::array<CF_Zero_Neumann, limite_num_parches>& g_zero_neumann
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
        std::vector<std::shared_ptr<Base>>& lista_parches_dinamicos,
        const std::array<CF_Dirichlet, limite_num_parches>& g_dirichlet,
        const std::array<CF_Zero_Neumann, limite_num_parches>& g_zero_neumann
    );


    // Funcion que retorna el tipo de CF. En caso de no encontrar el parche
    // retorna su nombre
    std::pair<std::string,int> que_tipo_es
    (
        const std::string& nombre,
        const std::array<CF_Dirichlet, limite_num_parches>& g_dirichlet,
        const std::array<CF_Zero_Neumann, limite_num_parches>& g_zero_neumann
    );

    // Funcion contenedor que asigna las condiciones de frontera para evitar
    // codigo muy gordo en el main
    void construir_condiciones_de_frontera
    (
        std::vector<Malla::Mallador::Parche>& Parches_norte,
        std::vector<Malla::Mallador::Parche>& Parches_sur,
        std::vector<Malla::Mallador::Parche>& Parches_este,
        std::vector<Malla::Mallador::Parche>& Parches_oeste,
        std::vector<double>& phi,
        int nx,
        std::vector<Dirichlet>& lista_parches_dirichlet,
        std::vector<std::shared_ptr<Base>>& lista_parches_dinamicos,
        const std::array<CF_Dirichlet, limite_num_parches>& g_dirichlet,
        const std::array<CF_Zero_Neumann, limite_num_parches>& g_zero_neumann
    );

    void construir_CF_flujo_de_masa
    (
        const Malla::Mallador& malla,
        std::vector<double>& u,
        std::vector<double>& v,
        std::vector<double>& velface,
        std::vector<double>& mdotstar
    );


} // Fin namespace Condicion_frontera

#endif //CONDICIONES_DE_FRONTERA_HPP
