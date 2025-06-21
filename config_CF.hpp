#pragma once
#include <string>
#include <array>

struct CF_Dirichlet {
    const std::string nombre;
    const double valor;
};

struct CF_Zero_Neumann {
    const std::string nombre;
    const std::string localizacion_fisica;
};

// Condiciones de frontera de Dirichlet
inline std::array<CF_Dirichlet, 3> g_dirichlet = {{
    {"norte", 100.0},
    {"sur", 55.0},
    {"este", 22.0}
    }
};

// Condiciones de frontera de Dirichlet
inline std::array<CF_Zero_Neumann, 4> g_zero_neumann = { {
        {"oeste", "oeste"}
    }
};
