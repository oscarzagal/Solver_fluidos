//
// Created by oscar on 8/06/25.
//

#include <stdexcept>

#include "condiciones_de_frontera.hpp"
#include "config_CF.hpp"

namespace Condicion_frontera {

    // Lista de inicializacion
    Dirichlet::Dirichlet
    (
        std::vector<int>& nodos_del_parche_,
        const double valor_,
        std::vector<double>& phi_
    ) :
    nodos_del_parche(nodos_del_parche_),
    valor(valor_),
    phi(phi_)
    {}

    void Dirichlet::aplicar() {

        for (const int& index : nodos_del_parche) {
            phi[index]=valor;
        }

    }

    // Lista de inicializacion
    Zero_Neumann::Zero_Neumann
    (
        std::vector<double>& phi_,
        const std::vector<int>& nodos_de_parche_,
        const std::string& frontera_fisica_,
        const int nx_
    ) :
    phi(phi_),
    nodos_del_parche(nodos_de_parche_),
    frontera_fisica(frontera_fisica_),
    nx(nx_)
    {}

    void Zero_Neumann::aplicar() {

        int desfase=0;

        if (frontera_fisica == "norte") {
            desfase=-nx;
        } else if (frontera_fisica == "sur") {
            desfase=nx;
        } else if (frontera_fisica == "este") {
            desfase=-1;
        } else if (frontera_fisica == "oeste") {
            desfase=1;
        } else {
            throw std::runtime_error("La frontera "
                + frontera_fisica + " no es un frontera");
        }


        for (const int& index : nodos_del_parche) {
            phi[index] = phi[index+desfase];
        }

    }

    std::shared_ptr<Base> Fabrica_de_CF::crear
    (
        const std::pair<std::string, int> &tipo,
        const Malla::Mallador::Parche &parche,
        std::vector<double> &phi,
        int nx,
        const std::array<CF_Zero_Neumann, limite_num_parches> &g_zero_neumann
    )
    {
        if (tipo.first == "zero_neumann") {

            return std::make_unique<Zero_Neumann>
            (
                phi,
                parche.obtener_nodos_del_parche,
                g_zero_neumann[tipo.second].localizacion_fisica,
                nx
            );

        }

        return {};

    }


    void asignar_condiciones_de_frontera
    (
        std::vector<Malla::Mallador::Parche>& parches,
        std::vector<double>& phi,
        const int& nx,
        std::vector<Dirichlet>& lista_dirichlet,
        std::vector<std::shared_ptr<Base>>& lista_parches_dinamicos,
        const std::array<CF_Dirichlet, limite_num_parches>& g_dirichlet,
        const std::array<CF_Zero_Neumann, limite_num_parches>& g_zero_neumann
    )
    {
        for (int i = 0; i < static_cast<int>(parches.size()); ++i) {
            std::pair<std::string,int> tipo = que_tipo_es
            (parches[i].obtener_nombre,
             g_dirichlet,
             g_zero_neumann);

            if (tipo.first == "dirichlet") {

                // NOTE: seria conveniente usar "emplace_back" para no crear
                // objetos temporales.
                Dirichlet parametros
                (
                    parches[i].obtener_nodos_del_parche,
                    g_dirichlet[tipo.second].valor,
                    phi
                );

                lista_dirichlet.push_back(parametros);

            } else if (tipo.first == "zero_neumann") {

                lista_parches_dinamicos.push_back(Fabrica_de_CF::crear(tipo,parches[i],phi,nx,g_zero_neumann));

            } else {
                throw std::runtime_error
                (
                    "El parche "
                    + tipo.first +
                    " no fue encontrado en config_CF.hpp"
                );
            }
        }
    }

    std::pair<std::string,int> que_tipo_es
    (
        const std::string& nombre,
        const std::array<CF_Dirichlet, limite_num_parches>& g_dirichlet,
        const std::array<CF_Zero_Neumann, limite_num_parches>& g_zero_neumann
    )
    {

        for (int i=0; i<static_cast<int>(g_dirichlet.size()); ++i) {
            if (g_dirichlet[i].nombre == nombre) {
                return {"dirichlet",i};
            }
        }

        for (int i=0; i<static_cast<int>(g_zero_neumann.size()); ++i) {
            if (g_zero_neumann[i].nombre == nombre) {
                return {"zero_neumann",i};
            }
        }

        return {nombre,{}};

    }

    void construir_condiciones_de_frontera
    (
        std::vector<Malla::Mallador::Parche>& Parches_norte,
        std::vector<Malla::Mallador::Parche>& Parches_sur,
        std::vector<Malla::Mallador::Parche>& Parches_este,
        std::vector<Malla::Mallador::Parche>& Parches_oeste,
        std::vector<double>& phi,
        const int nx,
        std::vector<Dirichlet>& lista_parches_dirichlet,
        std::vector<std::shared_ptr<Base>>& lista_parches_dinamicos,
        const std::array<CF_Dirichlet, limite_num_parches>& g_dirichlet,
        const std::array<CF_Zero_Neumann, limite_num_parches>& g_zero_neumann
    )
    {

        // NOTE: yo creo que con un "for" queda, o no???

        // Frontera norte
        asignar_condiciones_de_frontera
        (
            Parches_norte,
            phi,
            nx,
            lista_parches_dirichlet,
            lista_parches_dinamicos,
            g_dirichlet,
            g_zero_neumann
        );

        // Frontera sur
        asignar_condiciones_de_frontera
        (
            Parches_sur,
            phi,
            nx,
            lista_parches_dirichlet,
            lista_parches_dinamicos,
            g_dirichlet,
            g_zero_neumann
        );

        // Frontera este
        asignar_condiciones_de_frontera
        (
            Parches_este,
            phi,
            nx,
            lista_parches_dirichlet,
            lista_parches_dinamicos,
            g_dirichlet,
            g_zero_neumann
        );

        // Frontera oeste
        asignar_condiciones_de_frontera
        (
            Parches_oeste,
            phi,
            nx,
            lista_parches_dirichlet,
            lista_parches_dinamicos,
            g_dirichlet,
            g_zero_neumann
        );
    }



} // Fin namespace Condicion_frontera
