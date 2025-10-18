#ifndef MOLECULE_H
#define MOLECULE_H

// #define DEBUG_MODE

#include <string>
#include <vector>
#include <array>
#include <optional>
#include <cblas.h>
#include <omp.h>

#include "Atom.h"
#include "BasisFunction.h"
#include "Orbital.h"

class Molecule {
public:
    Molecule();
    Molecule(const std::string& molden_file_name);
    ~Molecule();

    static Molecule from_molden(const std::string& molden_file_name);

    void read_atoms_from_molden();
    void read_basis_functions_from_molden();
    void read_orbitals_from_molden();
    void update_molecule_limits();

    double scfp_density_at_point(const std::array<double, 3>& point) const;
    double orbital_value_at_point(const std::array<double, 3>& point, size_t number) const;

    bool has_density_matrix() const { return !density_matrix.empty(); }
    bool has_orbitals() const { return !orbitals.empty(); }

    // геттеры для границ молекулы
    double get_x_min() const { return x_min; }
    double get_x_max() const { return x_max; }
    double get_y_min() const { return y_min; }
    double get_y_max() const { return y_max; }
    double get_z_min() const { return z_min; }
    double get_z_max() const { return z_max; }

    std::vector<Atom> get_atoms() const { return atoms; }
    int get_number_of_orbitals() const { return orbitals.size(); }
    int get_number_of_occupied_orbitals() const { return number_of_occupied_orbitals; }
    const Orbital* get_orbitals() const { return orbitals.data(); }
    const BasisFunction* get_basis_functions() const { return basis_functions.data(); }


private:
    std::vector<Atom> atoms;
    std::vector<BasisFunction> basis_functions;
    std::vector<Orbital> orbitals;
    std::vector<std::vector<double>> density_matrix;

    std::string molden_data;

    double x_min, x_max, y_min, y_max, z_min, z_max;

    int number_of_occupied_orbitals;

    // Вспомогательные функции
    std::string read_file(const std::string& filename) const;
};

#endif // MOLECULE_H
