#include <iostream>
#include <vector>
#include <ranges>

int main() {

    std::vector<int> huevo = {1, 2, 3, 4, 5};
    std::vector<int> panza = {10, 20, 30, 40, 50};

    for (auto [i, valor, panza] : std::views::zip(std::views::iota(0), huevo, panza)) {

        std::cout << "i = " << i << " -> " << valor << " -> " << panza << "\n";

    }


    return 0;
}
