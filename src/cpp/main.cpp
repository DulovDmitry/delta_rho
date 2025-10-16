#include <iostream>
#include <chrono>

#include "Molecule.h"
#include "Grid.h"

int main() {
	std::string filename = "";
	auto mol = Molecule(filename);
	// auto mol = Molecule("../../H2.molden.input");

	auto p = mol.get_atoms()[0].get_position();

	auto grid = RegularOrthogonalGrid(mol);

	auto tStart = std::chrono::high_resolution_clock::now();

	// #pragma omp parallel for schedule(dynamic)
	// for (int i = 0; i < 1000; ++i){
	// 	double dens = mol.scfp_density_at_point(p);
	// }

	auto res = grid.calculate_scfp_values(mol);

	
	auto tEnd = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = tEnd - tStart;
	std::cout << "Time taken: " << elapsed.count() << "s\n";

	for (auto a : res) {
		std::cout << a << "\n";
	}


	tStart = std::chrono::high_resolution_clock::now();
	double dens = mol.scfp_density_at_point(p);
	tEnd = std::chrono::high_resolution_clock::now();
	elapsed = tEnd - tStart;
	std::cout << "Time taken: " << elapsed.count() << "s\n";
	std::cout << dens << "\n";


	return 0;
}