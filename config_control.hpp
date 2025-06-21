//
// Created by oscar on 13/06/25.
//

#pragma once

#include <string>

// Propiedades termofisicas
// Conductividad termica
inline constexpr double k = 20.0;

// Control
inline constexpr int num_iteraciones_max = 6000;

inline constexpr double tolerancia = 1e-10;

// Resolución del sistema de ecuaciones
// Solver lineal
inline std::string solver_elegido = "SOR";

// Factor de relajacion
inline constexpr double lambdaT = 1.5;
