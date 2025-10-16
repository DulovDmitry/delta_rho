#include "Grid.h"
#include <algorithm>
#include <numeric>
#include <iostream>
#include <stdexcept>

Grid::Grid() : scfp_values(){}

std::vector<double> Grid::calculate_scfp_values(const Molecule& molecule) {
    if (!molecule.has_density_matrix() && !molecule.has_orbitals()) {
        throw std::runtime_error("The molecule has neither density matrix nor orbitals");
    }

    scfp_values.clear();
    scfp_values.reserve(x.size());

    std::cout << "starting scfp calcilation\nNumber of grid points: " << x.size() << std::endl;

    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < x.size(); ++i) {
        std::array<double, 3> point = {x[i], y[i], z[i]};
        double val = molecule.scfp_density_at_point(point);
        scfp_values.push_back(val);
    }

    return scfp_values;
}

// -------- SphericalGrid --------

SphericalGrid::SphericalGrid(int phi_points, int theta_points, int r_points,
                             double radius, std::array<double, 3> center)
    : center_point(center),
      delta_phi(2 * M_PI / phi_points),
      delta_theta(M_PI / theta_points),
      delta_r(2 * radius / (2 * r_points - 1))
{
    double phi_initial = delta_phi / 2.0;
    double phi_final = 2 * M_PI - (delta_phi / 2.0);
    double theta_initial = delta_theta / 2.0;
    double theta_final = M_PI - (delta_theta / 2.0);
    double r_initial = delta_r / 2.0;

    std::vector<double> phi(phi_points), theta(theta_points), r(r_points);
    for (int i = 0; i < phi_points; ++i)
        phi[i] = phi_initial + i * (phi_final - phi_initial) / (phi_points - 1);
    for (int i = 0; i < theta_points; ++i)
        theta[i] = theta_initial + i * (theta_final - theta_initial) / (theta_points - 1);
    for (int i = 0; i < r_points; ++i)
        r[i] = r_initial + i * (radius - r_initial) / (r_points - 1);

    for (double rv : r) {
        for (double tv : theta) {
            for (double pv : phi) {
                double X = center[0] + rv * std::sin(tv) * std::cos(pv);
                double Y = center[1] + rv * std::sin(tv) * std::sin(pv);
                double Z = center[2] + rv * std::cos(tv);
                x.push_back(X);
                y.push_back(Y);
                z.push_back(Z);
                this->r.push_back(rv);
                this->theta.push_back(tv);
                this->phi.push_back(pv);
            }
        }
    }
}

// -------- RegularOrthogonalGrid --------

RegularOrthogonalGrid::RegularOrthogonalGrid(double x_length, double y_length, double z_length,
                                             int x_points_count, int y_points_count, int z_points_count,
                                             std::array<double, 3> center)
    : center_point(center)
{
    create_grid(x_length, y_length, z_length, x_points_count, y_points_count, z_points_count, center);

    for (auto p : x) {
        std::cout << p << "\n";
    }
}

RegularOrthogonalGrid::RegularOrthogonalGrid(const Molecule& molecule, int x_points_count, int y_points_count, int z_points_count)
{
    center_point = {
        (molecule.get_x_min() + molecule.get_x_max()) / 2.0,
        (molecule.get_y_min() + molecule.get_y_max()) / 2.0,
        (molecule.get_z_min() + molecule.get_z_max()) / 2.0
    };
    
    double x_margin = 2.0;
    double y_margin = 2.0;
    double z_margin = 2.0;

    double x_length = molecule.get_x_max() - molecule.get_x_min() + x_margin;
    double y_length = molecule.get_y_max() - molecule.get_y_min() + y_margin;
    double z_length = molecule.get_z_max() - molecule.get_z_min() + z_margin;

    create_grid(x_length, y_length, z_length, x_points_count, y_points_count, z_points_count, center_point);

    for (auto p : x) {
        std::cout << p << "\n";
    }
}

void RegularOrthogonalGrid::create_grid(double x_length, double y_length, double z_length,
                                        int x_points_count, int y_points_count, int z_points_count,
                                        std::array<double, 3> center)
{
    x_points.resize(x_points_count);
    y_points.resize(y_points_count);
    z_points.resize(z_points_count);

    for (int i = 0; i < x_points_count; ++i)
        x_points[i] = center[0] - x_length / 2.0 + i * x_length / (x_points_count - 1);
    for (int i = 0; i < y_points_count; ++i)
        y_points[i] = center[1] - y_length / 2.0 + i * y_length / (y_points_count - 1);
    for (int i = 0; i < z_points_count; ++i)
        z_points[i] = center[2] - z_length / 2.0 + i * z_length / (z_points_count - 1);

    for (double xv : x_points)
        for (double yv : y_points)
            for (double zv : z_points) {
                x.push_back(xv);
                y.push_back(yv);
                z.push_back(zv);
            }

    delta_x = x_length / (x_points_count - 1);
    delta_y = y_length / (y_points_count - 1);
    delta_z = z_length / (z_points_count - 1);
    delta_V = delta_x * delta_y * delta_z;
}

// Простая реализация численного интеграла по сумме (аналог Simpson)
double RegularOrthogonalGrid::scfp_integral() const {
    double sum = std::accumulate(scfp_values.begin(), scfp_values.end(), 0.0);
    return sum * delta_V;
}

// -------- IrregularOrthogonalGrid --------

IrregularOrthogonalGrid::IrregularOrthogonalGrid(double delta, double scaling_factor,
                                                 int points, std::array<double, 3> center)
    : center_point(center)
{
    std::vector<double> template_vals(points);
    template_vals[0] = delta;
    for (int i = 1; i < points; ++i)
        template_vals[i] = template_vals[i - 1] + std::pow(scaling_factor, i) * delta;

    std::vector<double> positive(template_vals.begin(), template_vals.end());
    std::vector<double> negative(template_vals.begin(), template_vals.end());
    std::reverse(negative.begin(), negative.end());
    for (double& v : negative) v *= -1;

    x_points.reserve(1 + positive.size() + negative.size());
    y_points.reserve(x_points.size());
    z_points.reserve(x_points.size());

    x_points.push_back(center[0]);
    y_points.push_back(center[1]);
    z_points.push_back(center[2]);

    for (double v : negative) {
        x_points.push_back(center[0] + v);
        y_points.push_back(center[1] + v);
        z_points.push_back(center[2] + v);
    }
    for (double v : positive) {
        x_points.push_back(center[0] + v);
        y_points.push_back(center[1] + v);
        z_points.push_back(center[2] + v);
    }

    std::sort(x_points.begin(), x_points.end());
    std::sort(y_points.begin(), y_points.end());
    std::sort(z_points.begin(), z_points.end());

    for (double xv : x_points)
        for (double yv : y_points)
            for (double zv : z_points) {
                x.push_back(xv);
                y.push_back(yv);
                z.push_back(zv);
            }
}

double IrregularOrthogonalGrid::scfp_integral() const {
    double sum = std::accumulate(scfp_values.begin(), scfp_values.end(), 0.0);
    return sum;
}
