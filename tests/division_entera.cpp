#include <iostream>
#include <vector>

int main() {

    constexpr int nx = 11;
    constexpr int ny = 11;

    // std::vector<int> vec = {0 ,1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    std::vector<int> vec = {110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120};

    std::vector<int> nodos_borrados;

    std::cout << "vec antes de borrar nodos esquina: \n";
    for (int i = 0; i < static_cast<int>(vec.size()); ++i) {
        std::cout << vec[i] << " ";
    }
    std::cout << "\n\n";

    for (int i = 0; i < static_cast<int>(vec.size()); ++i) {

        // El modulo obtiene el indice local "i"
        int index_i = vec[i] % nx;
        int index_j = vec[i] / nx;

        if ((index_i == 0      && index_j == 0)      ||
            (index_i == nx - 1 && index_j == ny - 1) ||
            (index_i == nx - 1 && index_j == 0)      ||
            (index_i == 0      && index_j == ny - 1)) {

            std::cout << "vec[" << i << "] = " << vec[i] << "\n";

            nodos_borrados.push_back(vec[i]);

            vec.erase(vec.begin() + i);

        }
    }


    std::cout << "vec despues de borrar nodos esquina: \n";
    for (int i = 0; i < static_cast<int>(vec.size()); ++i) {
        std::cout << vec[i] << " ";
    }
    std::cout << "\n\n";

    std::cout << "Nodos borrados: \n";
    for (int i = 0; i < static_cast<int>(nodos_borrados.size()); ++i) {
        std::cout << nodos_borrados[i] << " ";
    }


    return 0;
}
