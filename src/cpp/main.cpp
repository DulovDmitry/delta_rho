#include "Molecule.h"
#include "Grid.h"

int main() {
	auto mol = Molecule("../../H2.molden.input");
	// auto mol = Molecule("../../indene.molden.input");

	// std::cout << std::scientific << std::setprecision(12);
	// for (auto a : (mol.get_orbitals() + 1)->get_coefficients()) {
	// 	std::cout << a << "\n";
	// }

	auto p = mol.get_atoms()[0].get_position();
	auto grid = RegularOrthogonalGrid(mol, 100, 100, 100);

	auto tStart = std::chrono::high_resolution_clock::now();

	auto res = grid.calculate_scfp_values(mol);

	
	auto tEnd = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = tEnd - tStart;
	std::cout << "Time taken: " << elapsed.count() << "s\n";

	// grid.print_scfp_values();
	// grid.write_cube_file();

	std::cout << "Integral of electron density\nCalculated = " << grid.scfp_integral() << "\nExpected = " << mol.get_number_of_occupied_orbitals()*2 << std::endl;
	

	// tStart = std::chrono::high_resolution_clock::now();
	// double dens = mol.scfp_density_at_point(p);
	// tEnd = std::chrono::high_resolution_clock::now();
	// elapsed = tEnd - tStart;
	// std::cout << "Time taken: " << elapsed.count() << "s\n";
	// std::cout << dens << "\n";


	return 0;
}