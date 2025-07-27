//
// Created by oscar on 22/06/25.
//

#ifndef ECUACIONES_GOBERNANTES_HPP
#define ECUACIONES_GOBERNANTES_HPP

#include <vector>

#include "malla_por_bloques.hpp"

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

    class Momentum : public Base {
    public:
        // Constructor
        explicit Momentum(Malla::Mallador &);

        void ensamblar() override;

        A_coef obtener_coeficientes() override;

        A_coef A;

    private:

        void asignar_matriz(const A_coef& A_paso);

        Malla::Mallador& malla;

    };

} // Fin namespace Ecuaciones_gobernantes

#endif //ECUACIONES_GOBERNANTES_HPP
