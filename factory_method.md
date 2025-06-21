# Creación de la fabrica de condiciones de frontera

1. Clase base polimórfica:

```C++
// Condicion_frontera/Base.hpp
namespace Condicion_frontera {

class Base {
public:
    virtual void aplicar() = 0;
    virtual ~Base() {}
};

}
```

2. Clases derivadas:


```C++
// Dirichlet.hpp
class Dirichlet {
    std::vector<double>& phi_;
    std::vector<int> nodos_;
    double valor_;
public:
    Dirichlet(std::vector<double>& phi, const std::vector<int>& nodos, double valor);
    void aplicar(); // NO hace falta heredar de Base
};
```

```C++
// Zero_Neumann.hpp
#include "Base.hpp"

class Zero_Neumann : public Condicion_frontera::Base {
    std::vector<double>& phi_;
    std::vector<int> nodos_;
    std::string ubicacion_;
    int nx_;
public:
    Zero_Neumann(std::vector<double>& phi, const std::vector<int>& nodos, std::string ubicacion, int nx);
    void aplicar() override;
};
```

3. Fábrica solo para condiciones dinámicas

```C++
class FabricaCondiciones {
public:
    static std::unique_ptr<Condicion_frontera::Base> crear(
        const std::string& tipo,
        const Malla::Mallador::Parche& parche,
        std::vector<double>& phi,
        const std::map<std::string, std::string>& mapa_zero_neumann,
        int nx
    ) {
        if (tipo == "zero_neumann") {
            const auto& nombre = parche.obtener_nombre; // implemento el chatGPT, no yo
            auto it = mapa_zero_neumann.find(nombre);
            if (it == mapa_zero_neumann.end()) {
                                                 // parche.obtener_nombre
                throw std::out_of_range("Parche " + nombre + " no tiene configuración zero_neumann");
            }

            return std::make_unique<Condicion_frontera::Zero_Neumann>(
                phi,
                parche.obtener_nodos_del_parche,
                it->second,
                nx
            );
        }

        throw std::runtime_error("Tipo dinámico desconocido en fábrica: " + tipo);
    }
};
```

4. Función de asignación


```C++
std::vector<Condicion_frontera::Dirichlet> lista_estaticas;
std::vector<std::unique_ptr<Condicion_frontera::Base>> lista_dinamicas;

for (const auto& parche : parches) {
    const std::string& tipo = parche.obtener_tipo_de_condicion_de_frontera;
    const std::string& nombre = parche.obtener_nombre;

    if (tipo == "dirichlet") {
        auto it = Temp_frontera_para_dirichlet.find(nombre);
        if (it == Temp_frontera_para_dirichlet.end()) {
            throw std::out_of_range("Parche " + nombre + " no tiene configuración dirichlet");
        }

        lista_estaticas.emplace_back(phi, parche.obtener_nodos_del_parche, it->second);

    } else if (tipo == "zero_neumann") {
        lista_dinamicas.push_back(
            FabricaCondiciones::crear(tipo, parche, phi, ubicacion_fisica_zero_neumann, nx)
        );

    } else {
        throw std::runtime_error("Tipo de condición no reconocida: " + tipo);
    }
}
```

5. Aplicación

```C++
// Al inicio
for (auto& cf : lista_estaticas) {
    cf.aplicar();
}

// En cada paso de tiempo
for (auto& cf : lista_dinamicas) {
    cf->aplicar();
}
```
