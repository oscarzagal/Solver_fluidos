//
// Created by oscar on 13/06/25.
//

#pragma once

#include <string>

// Propiedades termofisicas
// Conductividad termica
inline constexpr double k = 20.0;

// Densidad
inline constexpr double rho = 1.0;

// Viscosidad cinemática
inline constexpr double nu = 1.0e-3;

// Control
inline constexpr int num_iteraciones_max = 6000;

inline constexpr double tolerancia = 1e-10;

// Resolución del sistema de ecuaciones
// Solver lineal
inline std::string solver_elegido_Temp = "SOR";
inline std::string solver_elegido_u = "SOR";
inline std::string solver_elegido_v = "SOR";
inline std::string solver_elegido_P = "SOR";

// Factores de relajacion
inline constexpr double lambda_T = 0.5;
inline constexpr double lambda_u = 0.5;
inline constexpr double lambda_v = 0.5;
inline constexpr double lambda_P = 0.5;
