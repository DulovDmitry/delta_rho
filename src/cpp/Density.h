#ifndef DENSITY_H
#define DENSITY_H

#include <vector>
#include <string>
#include <iostream>

#include "Grid.h"
#include "Molecule.h"

class Density
{
public:
	Density(Molecule* molecule, Grid* grid);
	~Density() = default;

	void calculate_values();
	void write_to_cube(std::string filename = "density.cub") const;
	double integral() const;

	std::vector<double> values() const { return values_; }

private:
	Molecule* molecule_;
	Grid* grid_;
	std::vector<double> values_;
};

#endif // DENSITY_H