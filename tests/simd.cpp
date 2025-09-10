// NOTE: PRUEBA FALLIDA PERO SE DEJA COMO REFERENCIA

#include <iostream>
#include <vector>
#include <experimental/simd>


int main() {
    using namespace std::experimental;

    std::vector<double> a = {1,2,3,4,5,6,7,8};
    std::vector<double> b = {10,20,30,40,50,60,70,80};
    std::vector<double> c;

    // simd<T> representa un "vector SIMD" del tipo más eficiente para T
    simd<double> va(a.data());
    simd<double> vb(b.data());

    simd<double> vc = va + vb; // suma vectorizada

    vc.copy_to(c.data(), vector_aligned);

    for (double x : c) {
        std::cout << x << " "; // imprime 11 22 33 44 55 66 77 88
    }
    std::cout << "\n";
}
