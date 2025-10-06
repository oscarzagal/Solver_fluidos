#include <iostream>
#include <fstream>
#include <sstream>
#include <cctype>
#include <vector>

int main() {

	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> phi;

	std::ifstream archivo_a_leer("T.dat");

	std::string linea;

	if (archivo_a_leer.is_open()) {
		while (getline(archivo_a_leer, linea)) {
			const char primer_caracter = linea[0];
			std::istringstream iss(linea);
			if (!std::isalpha(static_cast<unsigned char>(primer_caracter)) && !linea.empty()) {
				std::string x_string, y_string, phi_string;
				if (iss >> x_string >> y_string >> phi_string) {
					x.push_back(std::stod(x_string));
					y.push_back(std::stod(y_string));
					phi.push_back(std::stod(phi_string));
				}
			}

		}
		archivo_a_leer.close();
	} else {
		std::cerr << "No se pudo abrir el archivo, panzón \n";
	}

	for (int i = 0 ; i < static_cast<int>(x.size()) ; ++i) {
		std::cout << "x[" << i << "] = " << x[i] << "\n";
	}

	std::cout << "\n\n";

	for (int i = 0 ; i < static_cast<int>(y.size()) ; ++i) {
		std::cout << "y[" << i << "] = " << y[i] << "\n";
	}

	std::cout << "\n\n";

	for (int i = 0 ; i < static_cast<int>(phi.size()) ; ++i) {
		std::cout << "phi[" << i << "] = " << phi[i] << "\n";
	}



	return 0;
}
