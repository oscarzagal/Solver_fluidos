//
// Created by oscar on 22/06/25.
//

#ifndef ECUACIONES_GOBERNANTES_HPP
#define ECUACIONES_GOBERNANTES_HPP

#include "malla_por_bloques.hpp"
#include "variables_discretizacion.hpp"

#include <vector>

namespace Ecuaciones_gobernantes {

// Matriz A que almacena los coeficientes agrupados
struct A_coef {
    std::vector<double> ap, ae, aw, an, as, b;
};

class Base {
public:

    virtual void ensamblar() = 0;

    virtual A_coef obtener_coeficientes() = 0;

    virtual ~Base() = default;
};

class Energia : public Base {
public:

    // Constructor
    explicit Energia(Malla::Mallador &);

    void ensamblar() override;

    A_coef obtener_coeficientes() override;

    A_coef A;

private:

    void asignar_matriz(const A_coef& A_paso);

    Malla::Mallador& malla;

};

class Momentum : public Base{
public:
    // Constructor
    Momentum(double, Malla::Mallador &, std::vector<double> &);

    // Metodo que hace lo mismo que el metodo polimorfico "ensamblar"
    void unir_ecuacion();

    [[nodiscard]] A_coef obtener_coeficientes_para_u() const;
    [[nodiscard]] A_coef obtener_coeficientes_para_v() const;

    A_coef A_u, A_v;

    // Flujo de masa en las caras
    struct Mstar {
        // mstar_x(i+1,j) es para la cara derecha
        // mstar_x(i,j) es para la cara izquierda
        // mstar_y(i,j+1) es para la cara norte
        // mstar_y(i,j) es para la cara sur
        std::vector<double> mstar_x={}, mstar_y={};
    };

    // Nodos en la direccion "x" y "y"
    const int nx, ny;

    // Instancia de mstar
    Mstar mstar;

    // Viscosidad cinemática
    const double nu;

private:

    void asignar_matriz_para_u(const A_coef& A_paso);
    void asignar_matriz_para_v(const A_coef& A_paso);


    Malla::Mallador& malla;

    Malla::Mallador::Interpolacion inter;
    Grad_explicito gradP_explicito;
    fluxes_convectivos flux_conv;
    fluxes_difusivos flux_dif;

    // Campo de presion necesario para el armado de la ecuacion de momentum
    std::vector<double>& Pstar;

};

class Presion : public Base {
public:

    // Constructor
    explicit Presion(Malla::Mallador &);

    void ensamblar() override;

    A_coef obtener_coeficientes() override;

    A_coef A;

private:

    const int nx, ny;

};


} // Fin namespace Ecuaciones_gobernantes

#endif //ECUACIONES_GOBERNANTES_HPP
