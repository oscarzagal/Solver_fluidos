//
// Created by oscar on 10/09/25.
//

#include "variables_discretizacion.hpp"

// Implementacion de los constructores

Gradiente::Gradiente(const int nx, const int ny)
: grad_x_vol(nx * ny), grad_y_vol(nx * ny)
{}

fluxes_convectivos::fluxes_convectivos(const int nx, const int ny)
: fluxFConv_e(nx * ny), fluxFConv_w(nx * ny), fluxFConv_n(nx * ny), fluxFConv_s(nx * ny),
  fluxCConv_e(nx * ny), fluxCConv_w(nx * ny), fluxCConv_n(nx * ny), fluxCConv_s(nx * ny)
{}

fluxes_difusivos::fluxes_difusivos(const int nx, const int ny)
: fluxFDif_e(nx * ny), fluxFDif_w(nx * ny), fluxFDif_n(nx * ny), fluxFDif_s(nx * ny),
  fluxCDif_e(nx * ny), fluxCDif_w(nx * ny), fluxCDif_n(nx * ny), fluxCDif_s(nx * ny)
{}
