#ifndef GRID_H
#define GRID_H

#include <vector>
#include <array>
#include <iostream>
#include <string>
#include <cmath>
#include <stdexcept>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <numeric>

#include "Molecule.h"

class Grid {
public:
    Grid() = default;
    virtual ~Grid() = default;

    std::vector<double> x() const { return x_; }
    std::vector<double> y() const { return y_; }
    std::vector<double> z() const { return z_; }

protected:
    Molecule mol;

    std::vector<double> x_, y_, z_;
    std::array<double, 3> center_point_{};
};


class SphericalGrid : public Grid {
public:
    explicit SphericalGrid(int phi_points = 20, int theta_points = 10, int r_points = 10,
                           double radius = 3.0, std::array<double, 3> center = {0.0, 0.0, 0.0});

private:
    std::vector<double> r, phi, theta;
    double delta_phi, delta_theta, delta_r;
};


class OrthogonalGrid : public Grid {
public:
    OrthogonalGrid() = default;
    ~OrthogonalGrid() = default;

    std::vector<double> x_points() const { return x_points_; }
    std::vector<double> y_points() const { return y_points_; }
    std::vector<double> z_points() const { return z_points_; }
    int number_of_x_points() const { return number_of_x_points_; }
    int number_of_y_points() const { return number_of_y_points_; }
    int number_of_z_points() const { return number_of_z_points_; }


protected:
    std::vector<double> x_points_, y_points_, z_points_;
    int number_of_x_points_{}, number_of_y_points_{}, number_of_z_points_{};
};


class RegularOrthogonalGrid final : public OrthogonalGrid {
public:
    explicit RegularOrthogonalGrid(double x_length = 4.0, double y_length = 4.0, double z_length = 4.0,
                                   int x_points = 10, int y_points = 10, int z_points = 10,
                                   std::array<double, 3> center = {0.0, 0.0, 0.0});

    explicit RegularOrthogonalGrid(const Molecule& molecule, int x_points = 20, int y_points = 20, int z_points = 20);

    void create_grid(double x_length, double y_length, double z_length, int x_points_count, int y_points_count, int z_points_count, std::array<double, 3> center);

    double delta_x() const { return delta_x_; }
    double delta_y() const { return delta_y_; }
    double delta_z() const { return delta_z_; }
    double delta_V() const { return delta_V_; }

private:
    double delta_x_, delta_y_, delta_z_, delta_V_;
};


class IrregularOrthogonalGrid final : public OrthogonalGrid {
public:
    explicit IrregularOrthogonalGrid(double delta = 0.01, double scaling_factor = 1.3,
                                     int points = 20, std::array<double, 3> center = {0.0, 0.0, 0.0});

};

#endif // GRID_H
