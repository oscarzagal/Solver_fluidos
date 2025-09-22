std::ranges::transform:

```cpp
#include <ranges>
namespace rng = std::ranges;

auto make_pfm = [&](const Parche& src){
    Parches_Flujo_de_Masa pfm(nx, ny);
    pfm.obtener_nodos_del_parche = src.obtener_nodos_del_parche;
    pfm.obtener_nombre           = src.obtener_nombre;
    pfm.cortar_nodos_esquina();
    pfm.calcular_vector_normal_unitario();
    pfm.asignar_deltas(malla);
    pfm.calcular_desfase();
    pfm.tipo_de_CF = pfm.añadir_tipo_de_CF(g_dirichlet_u, g_zero_neumann_u);
    return pfm; // NRVO
};

auto build_vec = [&](auto& dst, const auto& src){
    dst.clear(); dst.reserve(src.size());
    rng::transform(src, std::back_inserter(dst), make_pfm);
};

build_vec(parches_norte_FM, Parches_norte);
build_vec(parches_sur_FM,   Parches_sur);
build_vec(parches_este_FM,  Parches_este);
build_vec(parches_oeste_FM, Parches_oeste);

```


