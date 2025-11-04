#include "Grid.h"

// -------- SphericalGrid --------

SphericalGrid::SphericalGrid(int phi_points, int theta_points, int r_points,
                             double radius, std::array<double, 3> center)
    : delta_phi(2 * M_PI / phi_points),
      delta_theta(M_PI / theta_points),
      delta_r(2 * radius / (2 * r_points - 1))
{
    center_point_ = center;

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
                x_.push_back(X);
                y_.push_back(Y);
                z_.push_back(Z);
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
                                                 : delta_x_(0), delta_y_(0), delta_z_(0), delta_V_(0)
{
    center_point_ = center;
    create_grid(x_length, y_length, z_length, x_points_count, y_points_count, z_points_count, center);
}

RegularOrthogonalGrid::RegularOrthogonalGrid(const Molecule& molecule,
                                             int x_points_count, int y_points_count, int z_points_count)
                                                 : delta_x_(0), delta_y_(0), delta_z_(0), delta_V_(0)
{
    center_point_ = {
        (molecule.get_x_min() + molecule.get_x_max()) / 2.0,
        (molecule.get_y_min() + molecule.get_y_max()) / 2.0,
        (molecule.get_z_min() + molecule.get_z_max()) / 2.0
    };
    
    double x_margin = 18.0;
    double y_margin = 18.0;
    double z_margin = 18.0;

    number_of_x_points_ = x_points_count;
    number_of_y_points_ = y_points_count;
    number_of_z_points_ = z_points_count;

    double x_length = molecule.get_x_max() - molecule.get_x_min() + x_margin;
    double y_length = molecule.get_y_max() - molecule.get_y_min() + y_margin;
    double z_length = molecule.get_z_max() - molecule.get_z_min() + z_margin;

    double x_points_per_bohr = x_points_count / x_length;
    double y_points_per_bohr = y_points_count / y_length;
    double z_points_per_bohr = z_points_count / z_length;

    std::cout << "Grid. Points per 1 Bohr (x, y, x): (" << x_points_per_bohr << ", " << y_points_per_bohr << ", " << z_points_per_bohr << ")\n";

    create_grid(x_length, y_length, z_length, x_points_count, y_points_count, z_points_count, center_point_);
}

void RegularOrthogonalGrid::create_grid(double x_length, double y_length, double z_length,
                                        int x_points_count, int y_points_count, int z_points_count,
                                        std::array<double, 3> center)
{
    x_points_.resize(x_points_count);
    y_points_.resize(y_points_count);
    z_points_.resize(z_points_count);

    delta_x_ = x_length / (x_points_count - 1);
    delta_y_ = y_length / (y_points_count - 1);
    delta_z_ = z_length / (z_points_count - 1);

    double x_min = center[0] - x_length / 2.0;
    double y_min = center[1] - y_length / 2.0;
    double z_min = center[2] - z_length / 2.0;

    for (int i = 0; i < x_points_count; ++i)
        x_points_[i] = x_min + i * delta_x_;
    for (int i = 0; i < y_points_count; ++i)
        y_points_[i] = y_min + i * delta_y_;
    for (int i = 0; i < z_points_count; ++i)
        z_points_[i] = z_min + i * delta_z_;

    for (double xv : x_points_)
        for (double yv : y_points_)
            for (double zv : z_points_) {
                x_.push_back(xv);
                y_.push_back(yv);
                z_.push_back(zv);
            }

    delta_V_ = delta_x_ * delta_y_ * delta_z_;
}


// -------- IrregularOrthogonalGrid --------

IrregularOrthogonalGrid::IrregularOrthogonalGrid(double delta, double scaling_factor,
                                                 int points, std::array<double, 3> center)
{
    center_point_ = center;

    std::vector<double> template_vals(points);
    template_vals[0] = delta;
    for (int i = 1; i < points; ++i)
        template_vals[i] = template_vals[i - 1] + std::pow(scaling_factor, i) * delta;

    std::vector<double> positive(template_vals.begin(), template_vals.end());
    std::vector<double> negative(template_vals.begin(), template_vals.end());
    std::reverse(negative.begin(), negative.end());
    for (double& v : negative) v *= -1;

    x_points_.reserve(1 + positive.size() + negative.size());
    y_points_.reserve(x_points_.size());
    z_points_.reserve(x_points_.size());

    x_points_.push_back(center[0]);
    y_points_.push_back(center[1]);
    z_points_.push_back(center[2]);

    for (double v : negative) {
        x_points_.push_back(center[0] + v);
        y_points_.push_back(center[1] + v);
        z_points_.push_back(center[2] + v);
    }
    for (double v : positive) {
        x_points_.push_back(center[0] + v);
        y_points_.push_back(center[1] + v);
        z_points_.push_back(center[2] + v);
    }

    std::sort(x_points_.begin(), x_points_.end());
    std::sort(y_points_.begin(), y_points_.end());
    std::sort(z_points_.begin(), z_points_.end());

    for (double xv : x_points_)
        for (double yv : y_points_)
            for (double zv : z_points_) {
                x_.push_back(xv);
                y_.push_back(yv);
                z_.push_back(zv);
            }
}