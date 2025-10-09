#include <iostream>
#include <vector>
#include <algorithm>

int main() {
    std::vector<double> numbers = {4.2, 7.1, 2.231, 8.90, 6.2, 1.1, 9.237, 3.234};

    // Find the iterator to the maximum element
    double max= *(std::max_element(numbers.begin(), numbers.end()));

    // Check if the vector is not empty and get the maximum value
    // if (max_it != numbers.end()) {
    std::cout << "The maximum element is: " << max << std::endl;
    // }

    return 0;
}   
