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
// NOTE: por el momento no se usa, pero se deja como referencia
// inline std::array<CF_Dirichlet, limite_num_parches> g_dirichlet_T = {{
//     {"norte_derecha", 0.0},
//     {"este_arriba", 0.0},
//     {"sur_izquierda", 100.0},
//     {"oeste_abajo", 100.0}
//     }
// };

// Parametros para la condicion de frontera tipo zero_neumann para T
// NOTE: por el momento no se usa, pero se deja como referencia
// inline std::array<CF_Zero_Neumann, limite_num_parches> g_zero_neumann_T = { {
//     {"oeste_arriba", "oeste"},
//     {"norte_izquierda", "norte"},
//     {"este_abajo", "este"},
//     {"sur_derecha","sur"}
//     }
// };

// NOTE: condiciones de frontera para el caso lid driven cavity

// Condiciones tipo dirichlet para la velocidad en u
inline std::array<CF_Dirichlet, limite_num_parches> g_dirichlet_u = {
    {
    {"norte", 1.0},
    {"sur", 0.0},
    {"este", 0.0},
    {"oeste", 0.0}
    }
};

// Condiciones tipo dirichlet para la velocidad en u
inline std::array<CF_Dirichlet, limite_num_parches> g_dirichlet_v = {
    {
    {"norte", 0.0},
    {"sur", 0.0},
    {"este", 0.0},
    {"oeste", 0.0}
    }
};

// Condiciones tipo zero_neumann para la velocidad en u
inline std::array<CF_Zero_Neumann, limite_num_parches> g_zero_neumann_u = {
    {

    }
};

// Condiciones tipo zero_neumann para la velocidad en v
inline std::array<CF_Zero_Neumann, limite_num_parches> g_zero_neumann_v = {
    {

    }
};


inline std::array<CF_Dirichlet,limite_num_parches> g_dirichlet_P = {{

}};

inline std::array<CF_Zero_Neumann,limite_num_parches> g_zero_neumann_P = {{
    {"norte", "norte"},
    {"sur", "sur"},
    {"este", "este"},
    {"oeste", "oeste"}
}};
