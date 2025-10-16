#ifndef GRID_H
#define GRID_H

#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <stdexcept>

#include "Molecule.h"

class Grid {
public:
    Grid();
    virtual ~Grid() = default;

    virtual double scfp_integral() const = 0;

    std::vector<double> scfp_values;
    std::vector<double> calculate_scfp_values(const Molecule& molecule);

protected:

    // координаты точек (для общих сеток)
    std::vector<double> x, y, z;
};

class SphericalGrid : public Grid {
public:
    SphericalGrid(int phi_points = 20, int theta_points = 10, int r_points = 10,
                  double radius = 3.0, std::array<double, 3> center = {0.0, 0.0, 0.0});

private:
    std::vector<double> r, phi, theta;
    std::array<double, 3> center_point;
    double delta_phi, delta_theta, delta_r;
};

class RegularOrthogonalGrid : public Grid {
public:
    RegularOrthogonalGrid(double x_length = 4.0, double y_length = 4.0, double z_length = 4.0,
                          int x_points = 10, int y_points = 10, int z_points = 10,
                          std::array<double, 3> center = {0.0, 0.0, 0.0});

    RegularOrthogonalGrid(const Molecule& molecule, int x_points = 20, int y_points = 20, int z_points = 20);

    void create_grid(double x_length, double y_length, double z_length, int x_points_count, int y_points_count, int z_points_count, std::array<double, 3> center);

    double scfp_integral() const override;

private:
    std::vector<double> x_points, y_points, z_points;
    double delta_x, delta_y, delta_z, delta_V;
    std::array<double, 3> center_point;
};

class IrregularOrthogonalGrid : public Grid {
public:
    IrregularOrthogonalGrid(double delta = 0.01, double scaling_factor = 1.3,
                            int points = 20, std::array<double, 3> center = {0.0, 0.0, 0.0});

    double scfp_integral() const override;

private:
    std::vector<double> x_points, y_points, z_points;
    std::array<double, 3> center_point;
};

#endif // GRID_H
