#include <iostream>
#include <chrono>

#include "Molecule.h"

int main() {
	std::string filename = "";
	auto mol = Molecule(filename);
	// auto mol = Molecule("../../H2.molden.input");

	auto p = mol.get_atoms()[0].get_position();

	// auto tStart = std::chrono::high_resolution_clock::now();
	// for (int i = 0; i < 1000; ++i){
	// 	double dens = mol.scfp_density_at_point(p);
	// }
	// auto tEnd = std::chrono::high_resolution_clock::now();
	// std::chrono::duration<double> elapsed = tEnd - tStart;
	// std::cout << "Time taken: " << elapsed.count() << "s\n";


	auto tStart = std::chrono::high_resolution_clock::now();
	double dens = mol.scfp_density_at_point(p);
	auto tEnd = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = tEnd - tStart;
	std::cout << "Time taken: " << elapsed.count() << "s\n";
	std::cout << dens << "\n";


	return 0;
}