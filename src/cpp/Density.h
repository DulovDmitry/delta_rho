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

	Density operator-(const Density& other) const;

	void write_to_cube(std::string filename = "density.cub") const;
	double integral() const;

	std::vector<double> alpha_density() const { return alpha_density_; }
	std::vector<double> beta_density() const { return beta_density_; }

private:
	void calculate_values(std::string spin);
	void calculate_alpha_density();
	void calculate_beta_density();

	Molecule* molecule_;
	Grid* grid_;
	std::vector<double> alpha_density_;
	std::vector<double> beta_density_;
	std::vector<double> electron_density_;
	std::vector<double> spin_density_;
};

#endif // DENSITY_H