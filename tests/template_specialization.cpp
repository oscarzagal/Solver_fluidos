#include <iostream>
#include <vector>

struct MiTipo {
    double x, y;
};

struct otroTipo {
    const std::vector<int>& nodos_parche;
    const double valor;
    std::vector<double>& phi;
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

template<>
class Procesar<otroTipo> {
public:
    void aplicar() {
        std::cout << "otroTipo: valor=" << otrotipo.valor << ", nodos_parche size=" << otrotipo.nodos_parche.size() << "\n";
        std::cout << "phi size: " << otrotipo.phi.size() << "\n";

        otrotipo.phi.clear();

        std::cout << "phi size luego del clear(): " << otrotipo.phi.size() << "\n";

        otrotipo.phi.push_back(3.14);
    }

    // Los parámetros del constructor deben coincidir con los del struct otroTipo
    Procesar<otroTipo>(const std::vector<int>& nodos_, const double val_, std::vector<double>& phi_)
    : otrotipo{nodos_, val_, phi_}
    {
        std::cout << "Constructor de Procesar<otroTipo>\n";
    }

    otroTipo otrotipo;

};


int main() {

    std::vector<int> nodos = {1, 2, 3};
    const double val = 5.5;
    std::vector<double> phi = {0.1, 0.2, 0.3};

    std::vector<int> nodos1 = {2, 3, 4, 7, 6};
    const double val1 = 8.5;
    std::vector<double> phi1 = {1.2, 1.2};

    std::cout << "main: phi1[0] antes del emplace_back(): " << phi1[0] << "\n";

    // NOTE: este seria como el "std::vector<Dirichlet>&"
    // o como el "std::vector<std::shared_ptr<Base>>&"
    // de la funcion "asignar_condiciones_de_frontera"
    // del archivo "condiciones_de_frontera.*"
    std::vector<Procesar<otroTipo>> otros;

    // Procesar<otroTipo> huevo(nodos,val,phi);
    // Procesar<otroTipo> chamba(nodos1,val1,phi1);

    // "push_back" requiere de un objeto, "emplace_back" no
    // otros.push_back(huevo);
    // otros.push_back(chamba);

    otros.emplace_back(nodos, val, phi);
    otros.emplace_back(nodos1, val1, phi1);

    otros[0].aplicar();
    otros[1].aplicar();

    std::cout << "main: phi[0] luego del push_back(): " << otros[1].otrotipo.phi[0] << "\n";
    std::cout << "main: phi1[0] luego del emplace_back(): " << phi1[0] << "\n";


    Procesar<int> p1;
    p1.operar(100);

    Procesar<MiTipo> p2;
    p2.operar(MiTipo{3.2, 7.8});
}
