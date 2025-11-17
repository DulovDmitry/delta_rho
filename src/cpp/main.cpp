#include "Molecule.h"
#include "Grid.h"
#include "Density.h"

int main() {
	// auto mol = Molecule("../../H2.molden.input");
	auto mol1 = Molecule("../../indene.molden.input");
	auto mol2 = Molecule("../../indene_radcat_vertical.molden.input");

	auto grid = RegularOrthogonalGrid(mol1, 150, 130, 100);

	const auto density1 = Density(&mol1, &grid);
	const auto density2 = Density(&mol2, &grid);

	auto density = density1 - density2;

	density1.write_to_cube("d1.cub");
	density2.write_to_cube("d2.cub");
	density.write_to_cube("d.cub");

	double integral = density.integral();

	std::cout << "Electron density integral = " << integral << std::endl;
	// std::cout << "Expected value = " << mol1.number_of_electrons() << std::endl;

	return 0;
}

void benchmark() {
	auto mol = Molecule("../../indene.molden.input");

	std::vector<int> grid_points;
	std::vector<double> integrals;
	std::vector<double> times;

	for (int N = 20; N <= 500; N += 20) {
		auto grid = RegularOrthogonalGrid(mol, N, N, N);

		auto tStart = std::chrono::high_resolution_clock::now();
		const auto density = Density(&mol, &grid);
		auto tEnd = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = tEnd - tStart;

		double val = density.integral();

		grid_points.push_back(N);
		integrals.push_back(val);
		times.push_back(elapsed.count());

		std::cout << "Grid: N = " << N << ", number of points = " << N*N*N << "\n";
		std::cout << "Integral value = " << val << "\n";
		std::cout << "Time taken: " << elapsed.count() << "s\n\n";
	}

	 for (size_t i = 0; i < integrals.size(); ++i) {
	    std::cout << grid_points[i] << " : " << integrals[i] << " : " << times[i] << "\n";
	}
}