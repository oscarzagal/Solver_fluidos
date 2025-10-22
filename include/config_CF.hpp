#pragma once
#include <string>
#include <array>

constexpr int limite_num_parches = 20;

struct CF_Dirichlet {
    const std::string nombre;
    const double valor;
};

struct CF_Zero_Neumann {
    const std::string nombre;
    const std::string localizacion_fisica; // norte, sur, este, oeste
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

// WARNING: LOS TIPOS DE PARCHE PARA "u" Y "v" DEBEN DE SER LOS MISMOS, los
// valores pueden cambiar.
// TODO: tirar una excepcion cuando lo de arriba no se cumpla
// NOTE: condiciones de frontera para el caso lid driven cavity

// Condiciones tipo dirichlet para la velocidad en u
inline std::array<CF_Dirichlet, limite_num_parches> g_dirichlet_u = {
    {
    {"norte", 0.0},
    {"oeste_abajo", 0.0},
    {"sur", 0.0},
    {"oeste_arriba", 1.5}
    }
};

// Condiciones tipo dirichlet para la velocidad en v
inline std::array<CF_Dirichlet, limite_num_parches> g_dirichlet_v = {
    {
    {"norte", 0.0},
    {"oeste_abajo", 0.0},
    {"sur", 0.0},
    {"oeste_arriba", 0.0}
    }
};

// Condiciones tipo zero_neumann para la velocidad en u
inline std::array<CF_Zero_Neumann, limite_num_parches> g_zero_neumann_u = {
    {
        {"este_abajo", "este"},
        {"este_arriba", "este"}
    }
};

// Condiciones tipo zero_neumann para la velocidad en v
inline std::array<CF_Zero_Neumann, limite_num_parches> g_zero_neumann_v = {
    {
        {"este_abajo", "este"},
        {"este_arriba", "este"}
    }
};


inline std::array<CF_Dirichlet,limite_num_parches> g_dirichlet_P = {{

    {"este_abajo", 0.0},
    {"este_arriba", 0.0},

}};

inline std::array<CF_Zero_Neumann,limite_num_parches> g_zero_neumann_P = {{
    {"norte", "norte"},
    {"sur", "sur"},
    {"oeste_abajo", "oeste"},
    {"oeste_arriba", "oeste"}
}};
