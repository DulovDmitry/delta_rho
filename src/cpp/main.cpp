#include "Molecule.h"
#include "Grid.h"

int main() {
	// auto mol = Molecule("../../H2.molden.input");
	auto mol = Molecule("../../indene.molden.input");

	// Benchmark
	// ******
	std::vector<int> grid_points;
	std::vector<double> integrals;
	std::vector<double> times;

	for (int N = 10; N <= 500; N += 10) {
		auto grid = RegularOrthogonalGrid(mol, N, N, N);

		auto tStart = std::chrono::high_resolution_clock::now();
		auto res = grid.calculate_scfp_values(mol);
		auto tEnd = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = tEnd - tStart;

		double val = grid.scfp_integral();

		grid_points.push_back(N);
		integrals.push_back(val);
		times.push_back(elapsed.count());

		std::cout << "Grid: N = " << N << ", number of points = " << N*N*N << "\n";
		std::cout << "Integral value = " << val << "\n";
		std::cout << "Time taken: " << elapsed.count() << "s\n\n";
	}

	return 0;
	// ******

	auto grid = RegularOrthogonalGrid(mol, 100, 100, 100);

	auto tStart = std::chrono::high_resolution_clock::now();
	auto res = grid.calculate_scfp_values(mol);
	auto tEnd = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = tEnd - tStart;
	std::cout << "Time taken: " << elapsed.count() << "s\n";

	// grid.print_scfp_values();
	// grid.write_cube_file();

	std::cout << "Integral of electron density\nCalculated = " << grid.scfp_integral() << "\nExpected = " << mol.get_number_of_occupied_orbitals()*2 << std::endl;
	

	// auto p = mol.get_atoms()[0].get_position();
	// tStart = std::chrono::high_resolution_clock::now();
	// double dens = mol.scfp_density_at_point(p);
	// tEnd = std::chrono::high_resolution_clock::now();
	// elapsed = tEnd - tStart;
	// std::cout << "Time taken: " << elapsed.count() << "s\n";
	// std::cout << dens << "\n";


	return 0;
}