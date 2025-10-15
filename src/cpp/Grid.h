#ifndef GRID_H
#define GRID_H

#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <stdexcept>

class Grid {
public:
    Grid();
    virtual ~Grid() = default;

    virtual double scfp_integral() const = 0;

    std::vector<double> scfp_values;
    std::vector<double> calculate_scfp_values(const Molecule& molecule);

protected:
    std::vector<double> scfp_values;

    // координаты точек (для общих сеток)
    std::vector<double> x, y, z;
};

class SphericalGrid : public Grid {
public:
    SphericalGrid(int phi_points = 20, int theta_points = 10, int r_points = 10,
                  double radius = 3.0, std::array<double, 3> center = {0.0, 0.0, 0.0});

private:
    std::vector<double> x, y, z;
    std::vector<double> r, phi, theta;
    std::array<double, 3> center_point;
    double delta_phi, delta_theta, delta_r;
};

class RegularOrthogonalGrid : public Grid {
public:
    RegularOrthogonalGrid(double x_length = 4.0, double y_length = 4.0, double z_length = 4.0,
                          int x_points = 10, int y_points = 10, int z_points = 10,
                          std::array<double, 3> center = {0.0, 0.0, 0.0});

    static RegularOrthogonalGrid from_molecule(double x_min, double x_max,
                                               double y_min, double y_max,
                                               double z_min, double z_max,
                                               int x_points = 20, int y_points = 20, int z_points = 20);

    double scfp_integral() const override;

private:
    std::vector<double> x_points, y_points, z_points;
    std::vector<double> x, y, z;
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
    std::vector<double> x, y, z;
    std::array<double, 3> center_point;
};

#endif // GRID_H
