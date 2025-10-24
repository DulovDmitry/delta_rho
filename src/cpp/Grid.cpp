#include "Grid.h"
#include <algorithm>
#include <numeric>
#include <iostream>
#include <stdexcept>


Grid::Grid() : scfp_values(), batch_size(0), number_of_batches(0) {}

std::vector<double> Grid::calculate_scfp_values(const Molecule& molecule) {
    // Old version

    // if (!molecule.has_density_matrix() && !molecule.has_orbitals()) {
    //     throw std::runtime_error("The molecule has neither density matrix nor orbitals");
    // }

    // scfp_values.clear();
    // scfp_values.reserve(x.size());

    // std::cout << "starting scfp calcilation\nNumber of grid points: " << x.size() << std::endl;

    // #pragma omp parallel for schedule(dynamic)
    // for (size_t i = 0; i < x.size(); ++i) {
    //     std::array<double, 3> point = {x[i], y[i], z[i]};
    //     double val = molecule.scfp_density_at_point(point);
    //     scfp_values.push_back(val);
    // }

    mol = molecule;

    const int M = molecule.get_number_of_orbitals();       // число орбиталей
    const int nocc = molecule.get_number_of_occupied_orbitals();    // число занятых орбиталей
    const int G = x.size();    // число точек в сетке

    // std::cout << "Number of molecular orbitals: " << M << std::endl;
    // std::cout << "Number of occupied molecular orbitals: " << nocc << std::endl;
    // std::cout << "Number of grid points: " << G << std::endl;

    // std::ofstream file("debug.txt");
    // file << "C_occ\n\n";

    // std::vector<double> C_occ; // size M * nocc
    // auto orbs = molecule.get_orbitals();
    // for (int i = 0; i < nocc; ++i) {
    //     // std::cout << "MO number " << i << ": " << (orbs + i)->get_occupancy() << "\n";
    //     auto coefs = (orbs + i)->get_coefficients();
    //     for (const auto &a : coefs){
    //         C_occ.push_back(a);
    //         // file << a << "\n";
    //     }
    // }

    std::vector<double> C_occ(M * nocc);

    auto orbs = molecule.get_orbitals();
    for (int mu = 0; mu < M; ++mu) {
        for (int i = 0; i < nocc; ++i) {
            const auto &coefs = (orbs + i)->get_coefficients(); // coefs[mu]
            C_occ[mu * nocc + i] = coefs[mu];
        }
    }

    // std::cout << "C_occ size = " << C_occ.size() << "\n";

    // Плотность (на всех точках)
    scfp_values.resize(G);

    batch_size = 5000;

    // std::cout << "batch_size = " << batch_size << "\n";

    // file << "\n\nChi\n\n";
    auto bf = molecule.get_basis_functions();

    // omp_set_num_threads(8);

    #pragma omp parallel for schedule(dynamic)
    for (int batch_start = 0; batch_start < G; batch_start += batch_size) {
        int B = std::min(batch_size, G - batch_start);

        std::vector<double> Chi(B * M); //size B * M
        std::vector<double> Psi(B * nocc, 0.0);

        // 1. Вычисляем значения базисных функций в точках батча

        for (int i = 0; i < B; ++i) {
            int point_index = batch_start + i;
            std::array<double, 3> point = {x[point_index], y[point_index], z[point_index]};

            for (int j = 0; j < M; ++j) {
                Chi[i * M + j] =(bf + j)->value_at_point(point);
            }
        }


        // std::cout << "Chi size = " << Chi.size() << "\n";
        // for (auto a : Chi) file << a << "\n";

        // 2. Умножаем Chi × C_occ → Psi
        //     Chi: (B×M)
        //     C_occ: (M×nocc)
        //     Psi: (B×nocc)
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    B, nocc, M,
                    1.0,
                    Chi.data(), M,
                    C_occ.data(), nocc,
                    0.0,
                    Psi.data(), nocc);

        // 3. Считаем плотность: rho = 2 * sum_i |Psi_i|^2
        for (int g = 0; g < B; ++g) {
            double val = 0.0;
            for (int i = 0; i < nocc; ++i) {
                double psi = Psi[g * nocc + i];
                val += psi * psi;
            }
            // std::cout << val << std::endl;
            scfp_values[batch_start + g] = 2.0 * val;
        }

        // Предположим:
        // Psi — вектор длины G * nocc (row-major: точки подряд, каждая строка — nocc орбиталей)
        // G — количество точек в сетке
        // nocc — число занятых орбиталей
        // dV — объем вокселя (в Bohr³)
    }

    // file.close();
    return scfp_values;
}

void Grid::print_scfp_values() {
    std::ofstream file("scfp_values.txt");
    for (size_t i = 0; i < x.size(); ++i) {
        file << x[i] << "\t" << y[i] << "\t" << z[i] << "\t" << scfp_values[i] << "\n";
    }
    file.close();
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
}

RegularOrthogonalGrid::RegularOrthogonalGrid(const Molecule& molecule, int x_points_count, int y_points_count, int z_points_count)
{
    center_point = {
        (molecule.get_x_min() + molecule.get_x_max()) / 2.0,
        (molecule.get_y_min() + molecule.get_y_max()) / 2.0,
        (molecule.get_z_min() + molecule.get_z_max()) / 2.0
    };
    
    double x_margin = 18.0;
    double y_margin = 18.0;
    double z_margin = 18.0;

    number_of_x_points = x_points_count;
    number_of_y_points = y_points_count;
    number_of_z_points = z_points_count;

    double x_length = molecule.get_x_max() - molecule.get_x_min() + x_margin;
    double y_length = molecule.get_y_max() - molecule.get_y_min() + y_margin;
    double z_length = molecule.get_z_max() - molecule.get_z_min() + z_margin;

    double x_points_per_bohr = x_points_count / x_length;
    double y_points_per_bohr = y_points_count / y_length;
    double z_points_per_bohr = z_points_count / z_length;

    std::cout << "Grid. Points per 1 Bohr (x, y, x): (" << x_points_per_bohr << ", " << y_points_per_bohr << ", " << z_points_per_bohr << ")\n";

    create_grid(x_length, y_length, z_length, x_points_count, y_points_count, z_points_count, center_point);
}

void RegularOrthogonalGrid::create_grid(double x_length, double y_length, double z_length,
                                        int x_points_count, int y_points_count, int z_points_count,
                                        std::array<double, 3> center)
{
    x_points.resize(x_points_count);
    y_points.resize(y_points_count);
    z_points.resize(z_points_count);

    delta_x = x_length / (x_points_count - 1);
    delta_y = y_length / (y_points_count - 1);
    delta_z = z_length / (z_points_count - 1);

    double x_min = center[0] - x_length / 2.0;
    double y_min = center[1] - y_length / 2.0;
    double z_min = center[2] - z_length / 2.0;

    for (int i = 0; i < x_points_count; ++i)
        x_points[i] = x_min + i * delta_x;
    for (int i = 0; i < y_points_count; ++i)
        y_points[i] = y_min + i * delta_y;
    for (int i = 0; i < z_points_count; ++i)
        z_points[i] = z_min + i * delta_z;

    for (double xv : x_points)
        for (double yv : y_points)
            for (double zv : z_points) {
                x.push_back(xv);
                y.push_back(yv);
                z.push_back(zv);
            }

    delta_V = delta_x * delta_y * delta_z;
}

// Простая реализация численного интеграла по сумме (аналог Simpson)
double RegularOrthogonalGrid::scfp_integral() const {
    double sum = std::accumulate(scfp_values.begin(), scfp_values.end(), 0.0);
    return sum * delta_V;
}

void RegularOrthogonalGrid::write_cube_file()
{
    std::ofstream cube("scfp_density.cub");
    cube << std::fixed << std::setprecision(6);

    if (!cube.is_open()) {
        throw std::runtime_error("Cannot open file");
    }

    // ───── Заголовок ─────
    cube << "CUBE FILE GENERATED BY SCF CALCULATION\n";
    cube << "Electron density\n";

    auto atoms = mol.get_atoms();

    // Найдём начало сетки и векторы шагов
    double x0 = x.front();
    double y0 = y.front();
    double z0 = z.front();

    // ───── Первая строка: число атомов и начало сетки ─────
    cube << std::setw(5) << atoms.size()
         << std::setw(12) << x0
         << std::setw(12) << y0
         << std::setw(12) << z0 << "\n";

    // ───── Следующие три строки: размерность и векторы шагов ─────
    cube << std::setw(5) << number_of_x_points
         << std::setw(12) << delta_x
         << std::setw(12) << 0.0
         << std::setw(12) << 0.0 << "\n";

    cube << std::setw(5) << number_of_y_points
         << std::setw(12) << 0.0
         << std::setw(12) << delta_y
         << std::setw(12) << 0.0 << "\n";

    cube << std::setw(5) << number_of_z_points
         << std::setw(12) << 0.0
         << std::setw(12) << 0.0
         << std::setw(12) << delta_z << "\n";

    // ───── Координаты атомов ─────
    for (const auto &a : atoms) {
        cube << std::setw(5) << a.get_charge()
             << std::setw(12) << 0.0     // формальный заряд (обычно 0)
             << std::setw(12) << a.get_position()[0]
             << std::setw(12) << a.get_position()[1]
             << std::setw(12) << a.get_position()[2] << "\n";
    }

    // ───── Значения плотности ─────
    // Порядок: x изменяется быстрее всего, затем y, затем z
    int idx = 0;
    for (int k = 0; k < number_of_z_points; ++k) {
        for (int j = 0; j < number_of_y_points; ++j) {
            for (int i = 0; i < number_of_x_points; ++i, ++idx) {
                cube << std::scientific << std::setprecision(5)
                     << std::setw(13) << scfp_values[idx];
                if ((idx + 1) % 6 == 0) cube << "\n"; // 6 значений на строку
            }
            // cube << "\n";
        }
    }

    cube.close();
    std::cout << "Cube file has been written" << std::endl;
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
