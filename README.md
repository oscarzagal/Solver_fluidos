# Solver de flujo de fluidos sencillos

Este solver posee la capacidad de solucionar problemas de flujo de fluidos
newtonianos, laminares, con transferencia de calor en geometrías rectangulares
con un número arbitrario de parches.

## jueves 21 agosto 2025

**NOTE**: el polimorfismo en tiempo de ejecución volvió horrible a mi código,
muy complicado y difícil de mantener, por lo que decidí reescribirlo usando
plantillas (polimorfismo en tiempo de compilación) y usando el paradigma
orientada a datos.

## lunes 8 agosto 2025

**NOTE**: Para cuando toque el desarrollo de la ecuación de corrección de presión:
- En Dirichlet se hace todo el calculo normal pero P_{b}^{'} (presión en la
frontera) es cero, haciendo que a_{N} sea cero del lado derecho pero no en su
contribución al coeficiente central a_{C}.
- En Zero Neumann se sabe cuál es el flujo de masa desde el principio, por lo
que no hay contribución al coeficiente central por parte de a_{N} porque este
nunca se forma.

## Template specialization

```cpp
struct MiTipo {
    double x, y;
};

template<typename T>
class Procesar {
public:
    void operar(const T& valor) {
        std::cout << "Genérico: " << valor << "\n";
    }
};


template<>
class Procesar<MiTipo> {
public:
    void operar(const MiTipo& obj) {
        std::cout << "MiTipo: x=" << obj.x << ", y=" << obj.y << "\n";
    }
};


int main() {
    Procesar<int> p1;
    p1.operar(100);

    Procesar<MiTipo> p2;
    p2.operar(MiTipo{3.2, 7.8});
}

```
