#include "Molecule.h"
#include "Grid.h"
#include "Density.h"

int main() {
	auto mol = Molecule("../../H2.molden.input");

	auto grid =IrregularOrthogonalGrid(0.01, 1.04, 50, {0,0,0});
	const auto density = Density(&mol, &grid);

	std::cout << "integral = " << density.integral() << std::endl;


	// проверка работы функции вычитания двух плотностей
	{
		auto mol1 = Molecule("../../indene.molden.input");
		auto mol2 = Molecule("../../indene_radcat_vertical.molden.input");

		auto p = mol1.get_atoms()[0].get_position();

		// auto grid = RegularOrthogonalGrid(10,10,10, 100, 100, 100, p);
		auto grid =IrregularOrthogonalGrid(0.01, 50, p, 20);

		const auto density1 = Density(&mol1, &grid);
		const auto density2 = Density(&mol2, &grid);

		auto density = density1 - density2;

		// density1.write_to_cube("d1.cub");
		// density2.write_to_cube("d2.cub");
		// density.write_to_cube("d.cub");
		density.transform();
		// density.write_to_cube("density_transformed_X.cub");

		std::cout << "density1 - density2 integral = " << density.integral() << "\n";
	}

	// попытка рассчитать вклад в электронную плотность от отдельных атомов
	/*{
		double sum = 0;
		for (int atom_number = 0; atom_number < mol1.get_atoms().size(); ++atom_number) {
			// const int atom_number = 1;
			auto atom_position = mol1.get_atoms()[atom_number].get_position();
			auto grid_for_atom = RegularOrthogonalGrid(10, 10, 10, 100, 100, 100, atom_position);

			const auto density1_atom = Density(&mol1, &grid_for_atom, atom_number);
			double integral_atom = density1_atom.integral();
			sum += integral_atom;

			std::cout << "Electron density of 0th atom integral = " << integral_atom << std::endl;
		}

		std::cout << "Full electron density integral = " << density1.integral() << "\n";
		// std::cout << "Electron density of atom = " << integral_atom << "\n";
		std::cout << "Sum of atomic density integrals = " << sum << "\n";
		std::cout << "Expected value = " << mol1.number_of_electrons() << std::endl;
	}*/

	return 0;
}

void benchmark_RegOgrid() {
	auto mol1 = Molecule("../../indene.molden.input");
	auto mol2 = Molecule("../../indene_radcat_vertical.molden.input");
	auto p = mol1.get_atoms()[0].get_position();

	std::vector<int> grid_points;
	std::vector<double> integrals;
	std::vector<double> times;

	for (int N = 20; N <= 500; N += 20) {
		auto grid =IrregularOrthogonalGrid(0.01, N, p, 20);

		auto tStart = std::chrono::high_resolution_clock::now();
		const auto density1 = Density(&mol1, &grid);
		const auto density2 = Density(&mol2, &grid);
		auto density = density1 - density2;
		double val = density.integral();
		auto tEnd = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = tEnd - tStart;


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

void benchmark_IrrOgrid() {
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