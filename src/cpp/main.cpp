#include "Molecule.h"
#include "Grid.h"
#include "Density.h"

int main() {
	// auto mol = Molecule("../../H2.molden.input");
	auto mol = Molecule("../../indene.molden.input");

	// Benchmark
	// ******
	// std::vector<int> grid_points;
	// std::vector<double> integrals;
	// std::vector<double> times;

	// for (int N = 20; N <= 500; N += 20) {
	// 	auto grid = RegularOrthogonalGrid(mol, N, N, N);

	// 	auto tStart = std::chrono::high_resolution_clock::now();
	// 	auto res = grid.calculate_scfp_values(mol);
	// 	auto tEnd = std::chrono::high_resolution_clock::now();
	// 	std::chrono::duration<double> elapsed = tEnd - tStart;

	// 	double val = grid.scfp_integral();

	// 	grid_points.push_back(N);
	// 	integrals.push_back(val);
	// 	times.push_back(elapsed.count());

	// 	std::cout << "Grid: N = " << N << ", number of points = " << N*N*N << "\n";
	// 	std::cout << "Integral value = " << val << "\n";
	// 	std::cout << "Time taken: " << elapsed.count() << "s\n\n";
	// }

	//  for (size_t i = 0; i < integrals.size(); ++i) {
    //     std::cout << grid_points[i] << " : " << integrals[i] << " : " << times[i] << "\n";
    // }

	// return 0;
	// ******

	auto grid = RegularOrthogonalGrid(mol, 100, 100, 100);

	const auto density = Density(&mol, &grid);
	density.write_to_cube();
	double integral = density.integral();

	std::cout << "Electron density integral = " << integral << std::endl;
	std::cout << "Expected value = " << mol.get_number_of_occupied_orbitals()*2 << std::endl;

	return 0;
}