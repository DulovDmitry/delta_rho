#include "Density.h"

Density::Density(Molecule* molecule, Grid* grid)
	: molecule_(molecule),
	  grid_(grid)
{
    const size_t G = grid_->x().size();
    alpha_density_.resize(G);
    beta_density_.resize(G);
    spin_density_.resize(G);
    electron_density_.resize(G);

    if (molecule_->number_of_occupied_alpha_orbitals()) {
        calculate_alpha_density();
    }
    if (molecule_->number_of_occupied_beta_orbitals()) {
        calculate_beta_density();
    }

    if (molecule_->number_of_occupied_alpha_orbitals() && molecule_->number_of_beta_orbitals()) {
        for (size_t i = 0; i < G; ++i) {
            spin_density_[i] = alpha_density_[i] - beta_density_[i];
            electron_density_[i] = alpha_density_[i] + beta_density_[i];
        }
    } else {
        for (size_t i = 0; i < G; ++i) {
            electron_density_[i] = alpha_density_[i];
        }
    }
}

void Density::calculate_alpha_density() {
    std::cout << "Alpha density values calculation" << std::endl;
    auto tStart = std::chrono::high_resolution_clock::now();
    calculate_values("Alpha");
    auto tEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = tEnd - tStart;
    std::cout << "Time taken: " << elapsed.count() << "s\n";
}

void Density::calculate_beta_density() {
    std::cout << "Beta density values calculation" << std::endl;
    auto tStart = std::chrono::high_resolution_clock::now();
    calculate_values("Beta");
    auto tEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = tEnd - tStart;
    std::cout << "Time taken: " << elapsed.count() << "s\n";
}

void Density::calculate_values(std::string spin) {
    const int M = (spin == "Alpha" ? molecule_->number_of_alpha_orbitals() : molecule_->number_of_beta_orbitals());       // number of orbitals
    const int nocc = (spin == "Alpha" ? molecule_->number_of_occupied_alpha_orbitals() : molecule_->number_of_occupied_beta_orbitals());    // number of occupied orbitals
    const size_t G = grid_->x().size();    // number of grid points

    std::vector<double>& values = (spin == "Alpha" ? alpha_density_ : beta_density_);

    // std::cout << "Number of molecular orbitals: " << M << std::endl;
    // std::cout << "Number of occupied molecular orbitals: " << nocc << std::endl;
    // std::cout << "Number of grid points: " << G << std::endl;

    const auto bf = molecule_->get_basis_functions();
    const auto orbs = (spin == "Alpha" ? molecule_->get_alpha_orbitals() : molecule_->get_beta_orbitals());

    std::vector<double> C_occ(M * nocc);
    for (int mu = 0; mu < M; ++mu) {
        for (int i = 0; i < nocc; ++i) {
            const auto &coefs = (orbs + i)->get_coefficients(); // coefs[mu]
            C_occ[mu * nocc + i] = coefs[mu];
        }
    }

    std::vector<double> occupancies;
    occupancies.resize(nocc);
    for (int i = 0; i < nocc; ++i) {
        occupancies[i] = (orbs+i)->get_occupancy();
        // std::cout << i << ": " << occupancies[i] << std::endl;
    }

    // omp_set_num_threads(8);

    const std::vector<double> x = grid_->x();
    const std::vector<double> y = grid_->y();
    const std::vector<double> z = grid_->z();

    const size_t batch_size = 5000;

    #pragma omp parallel for schedule(dynamic)
    for (int batch_start = 0; batch_start < G; batch_start += batch_size) {
        size_t B = std::min(batch_size, G - batch_start);

        std::vector<double> Chi(B * M); //size B * M
        std::vector<double> Psi(B * nocc, 0.0);

        // 1. Calculation of basis functions values on grid points
        for (size_t i = 0; i < B; ++i) {
            int point_index = batch_start + i;
            std::array<double, 3> point = {x[point_index], y[point_index], z[point_index]};

            for (int j = 0; j < M; ++j) {
                Chi[i * M + j] =(bf + j)->value_at_point(point);
            }
        }

        // 2. Multiplying Chi × C_occ → Psi
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

        // 3. Calculation of density: rho = 2 * sum_i |Psi_i|^2
        for (size_t g = 0; g < B; ++g) {
            double val = 0.0;
            for (int i = 0; i < nocc; ++i) {
                double psi = Psi[g * nocc + i];
                val += occupancies[i] * psi * psi;
            }
            values[batch_start + g] = val;
        }
    }
}

void Density::write_to_cube(std::string filename) const
{
    // Check if the grid_ is RegularOrthogonalGrid
    const auto* rogrid = dynamic_cast<RegularOrthogonalGrid*>(grid_);
    if (rogrid == nullptr) {
		throw std::runtime_error("Class Density. Writing to a cube file is available only for density based on regular othogonal grid");
	}

    if (filename.find(".cub") == std::string::npos) {
		filename += ".cub";
	}

    std::ofstream cube(filename);
    cube << std::fixed << std::setprecision(6);

    if (!cube.is_open()) {
        throw std::runtime_error("Cannot open file");
    }

    // ───── Заголовок ─────
    cube << "CUBE FILE GENERATED BY DELTA_RHO PROGRAM\n";
    cube << "Electron density\n";

    const auto atoms = molecule_->get_atoms();

    // Найдём начало сетки и векторы шагов
    const double x0 = grid_->x()[0];
    const double y0 = grid_->y()[0];
    const double z0 = grid_->z()[0];

    // ───── Первая строка: число атомов и начало сетки ─────
    cube << std::setw(5) << atoms.size()
         << std::setw(12) << x0
         << std::setw(12) << y0
         << std::setw(12) << z0 << "\n";

    // ───── Следующие три строки: размерность и векторы шагов ─────
    cube << std::setw(5) << rogrid->number_of_x_points()
         << std::setw(12) << rogrid->delta_x()
         << std::setw(12) << 0.0
         << std::setw(12) << 0.0 << "\n";

    cube << std::setw(5) << rogrid->number_of_y_points()
         << std::setw(12) << 0.0
         << std::setw(12) << rogrid->delta_y()
         << std::setw(12) << 0.0 << "\n";

    cube << std::setw(5) << rogrid->number_of_z_points()
         << std::setw(12) << 0.0
         << std::setw(12) << 0.0
         << std::setw(12) << rogrid->delta_z() << "\n";

    // ───── Координаты атомов ─────
    for (const auto &a : atoms) {
        cube << std::setw(5) << a.get_charge()
             << std::setw(12) << 0.0     // формальный заряд
             << std::setw(12) << a.get_position()[0]
             << std::setw(12) << a.get_position()[1]
             << std::setw(12) << a.get_position()[2] << "\n";
    }

    // ───── Значения плотности ─────
    // Порядок: x изменяется быстрее всего, затем y, затем z
    int idx = 0;
    for (int k = 0; k < rogrid->number_of_z_points(); ++k) {
        for (int j = 0; j < rogrid->number_of_y_points(); ++j) {
            for (int i = 0; i < rogrid->number_of_x_points(); ++i, ++idx) {
                cube << std::scientific << std::setprecision(5)
                     << std::setw(13) << electron_density_[idx];
                if ((idx + 1) % 6 == 0) cube << "\n"; // 6 значений на строку
            }
            // cube << "\n";
        }
    }

    cube.close();
    std::cout << "Cube file has been written" << std::endl;
}

double Density::integral() const {
    if ( dynamic_cast<RegularOrthogonalGrid*>(grid_) != nullptr ) {
        const auto rogrid = dynamic_cast<RegularOrthogonalGrid*>(grid_);
        double sum = std::accumulate(electron_density_.begin(), electron_density_.end(), 0.0);
        return sum * rogrid->delta_V();
    }
    else {
        return 0.0;
    }
}

Density Density::operator-(const Density& other) const
{
    // Check grids
    if (this->grid_ != other.grid_) {
        throw std::runtime_error("Error: Cannot subtract density objects calculated on different grids");
    }

    Density result(*this);

    const size_t G = result.electron_density_.size();
    if (other.electron_density_.size() != G) {
        throw std::runtime_error("Error: Density grids have different sizes");
    }

    for (size_t i = 0; i < G; ++i) {
        result.alpha_density_[i] = this->alpha_density_[i] - other.alpha_density_[i];
        result.beta_density_[i] = this->beta_density_[i] - other.beta_density_[i];
        result.electron_density_[i] = this->electron_density_[i] - other.electron_density_[i];
        result.spin_density_[i] = this->spin_density_[i] - other.spin_density_[i];
    }

    return result;
}

