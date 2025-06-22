#pragma once
#include <string>
#include <array>

constexpr int limite_num_parches = 10;

struct CF_Dirichlet {
    const std::string nombre;
    const double valor;
};

struct CF_Zero_Neumann {
    const std::string nombre;
    const std::string localizacion_fisica;
};

// Parametros para la condicion de frontera tipo Dirichlet para T
inline std::array<CF_Dirichlet, limite_num_parches> g_dirichlet_T = {{
    {"norte", 100.0},
    {"sur", 0.0},
    {"este", 0.0},
    {"oeste",0.0}
    }
};

// Parametros para la condicion de frontera tipo zero_neumann para T
inline std::array<CF_Zero_Neumann, limite_num_parches> g_zero_neumann_T = { {
    }
};
