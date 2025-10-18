#include <iostream>
#include <chrono>

#include "Molecule.h"
#include "Grid.h"

int main() {
	std::string filename = "/mnt/d/lambda_project/BSF-33_opt_b3lyp.molden.input";
	auto mol = Molecule(filename);
	// auto mol = Molecule("../../H2.molden.input");

	for (auto a : (mol.get_orbitals() + 1)->get_coefficients()) {
		std::cout << a << "\n";
	}

	auto p = mol.get_atoms()[0].get_position();

	auto grid = RegularOrthogonalGrid(mol, 50, 50, 50);

	auto tStart = std::chrono::high_resolution_clock::now();

	// #pragma omp parallel for schedule(dynamic)
	// for (int i = 0; i < 1000; ++i){
	// 	double dens = mol.scfp_density_at_point(p);
	// }

	auto res = grid.calculate_scfp_values(mol);

	
	auto tEnd = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = tEnd - tStart;
	std::cout << "Time taken: " << elapsed.count() << "s\n";

	grid.print_scfp_values();
	grid.write_cube_file();
	
	// for (auto a : res) {
	// 	std::cout << a << "\n";
	// }


	tStart = std::chrono::high_resolution_clock::now();
	double dens = mol.scfp_density_at_point(p);
	tEnd = std::chrono::high_resolution_clock::now();
	elapsed = tEnd - tStart;
	std::cout << "Time taken: " << elapsed.count() << "s\n";
	std::cout << dens << "\n";


	return 0;
}