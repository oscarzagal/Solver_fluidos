//
// Created by oscar on 13/06/25.
//

#pragma once

#include <string>

// VARIABLE GLOBAL: Número de campos a resolver
inline constexpr int NUM_CAMPOS = 5;

// Propiedades termofisicas

// VARIABLE GLOBAL: Densidad
inline constexpr double rho = 1.0;

// VARIABLE GLOBAL: Viscosidad cinemática
inline constexpr double nu = 2.5e-3;

// Control

// VARIABLE GLOBAL: Paso de tiempo
inline constexpr double delta_t = 0.5;

// VARIABLE GLOBAL: Numero maximo de iteraciones
inline constexpr int num_iteraciones_max = 600000;

// VARIABLE GLOBAL: Tolerancia
inline constexpr double tolerancia = 1e-11;

// Resolución del sistema de ecuaciones

// VARIABLE GLOBAL: Solver lineal para u
inline std::string solver_elegido_u = "SOR";

// VARIABLE GLOBAL: Solver lineal para v
inline std::string solver_elegido_v = "SOR";

// VARIABLE GLOBAL: Solver lineal para P
inline std::string solver_elegido_P = "SOR";

// Factores de relajacion para las ecuaciones de gobierno

// VARIABLE GLOBAL: Factor de relajacion para U
inline constexpr double lambda_Vel = 0.6;

// VARIABLE GLOBAL: Factor de relajacion para P
inline constexpr double lambda_P = 0.3;

// Factores de relajacion para las ecuaciones de gobierno para el solver lineal

// VARIABLE GLOBAL: Factor de relajacion para u usado en solver lineal
inline constexpr double lambda_Vel_SL = 1.0;

// VARIABLE GLOBAL: Factor de relajacion para P usado en solver lineal
// (en este caso aun no se que valor para phi_old se puede utilizar como
// parametro para el de la iteracion anterior, por lo que configuro a esta
// variable como 1.0 para ignorar el efecto de phi_old)
inline constexpr double lambda_P_SL = 1.0;
