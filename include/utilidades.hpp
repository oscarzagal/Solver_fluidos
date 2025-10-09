//
// Created by oscar on 26/07/25.
//

#ifndef UTILIDADES_HPP
#define UTILIDADES_HPP

// Funcion para calcular los indices absolutos, por ejemplo, los indices para
// las celdas interiores
inline int idx(const int i, const int j, const int nx) {
    return i + nx * j;
}

// Calculo de valores interpolados
inline double interpolar(const double psiC, const double psiN, const double gN) {
    return gN * psiN + (1.0 - gN) * psiC;
}

#endif //UTILIDADES_HPP
