#include "Molecule.h"

#include <fstream>
#include <sstream>
#include <regex>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <iterator>

Molecule::Molecule()
    : x_min(0), x_max(0), y_min(0), y_max(0), z_min(0), z_max(0), number_of_occupied_orbitals(0) {}

Molecule::Molecule(const std::string& molden_file_name) {
    std::cout << "** Loading a molden file ..." << std::endl;
    molden_data = read_file(molden_file_name);
    
    std::cout << "** Reading atom data ..." << std::endl;
    read_atoms_from_molden();
    #ifdef DEBUG_MODE
    for (auto a : atoms) {
        std::cout << a.repr() << "\n\n";
    }
    std::cout << std::endl;
    #endif

    std::cout << "** Reading GTO data ..." << std::endl;
    read_basis_functions_from_molden();
    #ifdef DEBUG_MODE
    for (auto b : basis_functions) {
        std::cout << b.repr() << "\n\n";
    }
    std::cout << std::endl;
    
    #endif

    std::cout << "** Reading MO data ..." << std::endl;
    read_orbitals_from_molden();
    #ifdef DEBUG_MODE
    for (auto o : orbitals) {
        std::cout << o.repr() << std::endl;
    }
    std::cout << std::endl;
    #endif

    std::cout << "A molden file has been read!" << std::endl;
}

Molecule::~Molecule() {}

std::string Molecule::read_file(const std::string& filename) const {
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Cannot open molden file: " + filename);

    std::ostringstream ss;
    ss << file.rdbuf();
    return ss.str();
}

void Molecule::read_atoms_from_molden() {
    const std::string start_block = "[Atoms] AU";
    const std::string end_block = "[GTO]";
    size_t start_index = molden_data.find(start_block);
    if (start_index == std::string::npos)
        throw std::runtime_error("Class Molecule. There is no [Atoms] block in the file");

    start_index += start_block.size();
    size_t end_index = molden_data.find(end_block, start_index);
    std::string fragment = (end_index == std::string::npos)
        ? molden_data.substr(start_index)
        : molden_data.substr(start_index, end_index - start_index);

    std::istringstream lines(fragment);
    std::string line;
    while (std::getline(lines, line)) {
        if (line.empty()) continue;

        std::istringstream ls(line);
        std::string label;
        int number, charge;
        double x, y, z;
        std::vector<std::string> tokens;
        std::string token;
        while (ls >> token) tokens.push_back(token);
        if (tokens.size() != 6)
            throw std::runtime_error("Class Molecule. Unrecognized atom line in molden file");

        label = tokens[0];
        number = std::stoi(tokens[1]);
        charge = std::stoi(tokens[2]);
        x = std::stod(tokens[3]) * 0.529177249;
        y = std::stod(tokens[4]) * 0.529177249;
        z = std::stod(tokens[5]) * 0.529177249;

        atoms.emplace_back(std::array<double, 3>{x, y, z}, label, number, charge);
    }

    update_molecule_limits();
}

void Molecule::update_molecule_limits() {
    if (atoms.empty()) return;

    auto first = atoms.front().get_position();
    x_min = x_max = first[0];
    y_min = y_max = first[1];
    z_min = z_max = first[2];

    for (const auto& atom : atoms) {
        const auto& pos = atom.get_position();
        x_min = std::min(x_min, pos[0]);
        x_max = std::max(x_max, pos[0]);
        y_min = std::min(y_min, pos[1]);
        y_max = std::max(y_max, pos[1]);
        z_min = std::min(z_min, pos[2]);
        z_max = std::max(z_max, pos[2]);
    }
}

void Molecule::read_basis_functions_from_molden() {
    const std::string start_block = "[GTO]";
    const std::string end_block = "[5D]";
    size_t start_index = molden_data.find(start_block);
    if (start_index == std::string::npos)
        throw std::runtime_error("Class Molecule. There is no [GTO] block in the file");

    start_index += start_block.size();
    size_t end_index = molden_data.find(end_block, start_index);
    std::string fragment = (end_index == std::string::npos)
        ? molden_data.substr(start_index)
        : molden_data.substr(start_index, end_index - start_index);

    std::regex block_split(R"(\n\s*\n)");
    std::sregex_token_iterator it(fragment.begin(), fragment.end(), block_split, -1);
    std::sregex_token_iterator end;

    std::vector<std::string> blocks(it, end);

    for (size_t block_index = 0; block_index < blocks.size(); ++block_index) {
        const auto& block = blocks[block_index];
        std::regex gto_split(R"((?=[spdf]))");
        std::sregex_token_iterator git(block.begin(), block.end(), gto_split, -1);
        std::vector<std::string> gto_blocks(git, end);

        for (size_t j = 1; j < gto_blocks.size(); ++j) {
            std::istringstream GTO_stream(gto_blocks[j]);
            std::string line;
            std::vector<std::string> lines;
            while (std::getline(GTO_stream, line)) {
                if (!line.empty()) lines.push_back(line);
            }

            if (lines.empty()) continue;
            std::istringstream header(lines[0]);
            std::string shell;
            header >> shell;

            std::vector<double> exponents, coefficients;
            for (size_t i = 1; i < lines.size(); ++i) {
                std::istringstream val_line(lines[i]);
                double exp, coef;
                if (val_line >> exp >> coef) {
                    exponents.push_back(exp);
                    coefficients.push_back(coef);
                }
            }

            auto position = atoms[block_index].get_position();
            if (shell == "s") {
                basis_functions.emplace_back(coefficients, exponents, shell, position);
            } else if (shell == "p") {
                basis_functions.emplace_back(coefficients, exponents, shell, position, "x");
                basis_functions.emplace_back(coefficients, exponents, shell, position, "y");
                basis_functions.emplace_back(coefficients, exponents, shell, position, "z");
            } else if (shell == "d") {
                basis_functions.emplace_back(coefficients, exponents, shell, position, "z2");
                basis_functions.emplace_back(coefficients, exponents, shell, position, "xz");
                basis_functions.emplace_back(coefficients, exponents, shell, position, "yz");
                basis_functions.emplace_back(coefficients, exponents, shell, position, "x2y2");
                basis_functions.emplace_back(coefficients, exponents, shell, position, "xy");
            }
        }
    }
}

void Molecule::read_orbitals_from_molden() {
    std::cout << "begin\n";
    const std::string start_block = "[MO]";

    std::cout << "finding MO block\n";
    size_t start_index = molden_data.find(start_block);
    if (start_index == std::string::npos)
        throw std::runtime_error("Class Molecule. There is no [MO] block in the file");

    std::cout << "creating substring\n";
    start_index += start_block.size();
    std::string fragment = molden_data.substr(start_index);
    // std::cout << fragment << std::endl;

    std::cout << "dividing substring into blocks\n";
    // std::regex mo_split(R"((?=Sym=))");
    // std::sregex_token_iterator it(fragment.begin(), fragment.end(), mo_split, -1);
    // std::sregex_token_iterator end;
    // std::vector<std::string> blocks(it, end);
    // for (auto s : blocks) std::cout << s << std::endl;

    std::string delimiter = "Sym=";
    std::vector<std::string> blocks;

    size_t start = 0;
    size_t end = 0;

    while (end != std::string::npos) {
        end = fragment.find(delimiter, start);
        std::string block = delimiter + fragment.substr(start, end - start);
        // std::cout << block << std::endl;

        blocks.push_back(block);
        start = end + delimiter.length();
    }

    std::cout << "Size of std::vector<std::string> blocks = " << blocks.size() << std::endl;

    std::cout << "reading blocks\n";
    int n = 0;
    for (size_t i = 1; i < blocks.size(); ++i) {
        std::istringstream block(blocks[i]);
        std::string line;
        std::vector<std::string> lines;
        while (std::getline(block, line)) {
            if (!line.empty()) lines.push_back(line);
        }
        if (lines.size() < 5) continue;

        double occ = 0.0;
        {
            std::string occ_line = lines[3];
            size_t pos = occ_line.find("=");
            if (pos != std::string::npos)
                occ = std::stod(occ_line.substr(pos + 1));
        }

        if (occ != 0) {
            n++;
        }

        std::vector<double> coefs;
        for (size_t j = 4; j < lines.size(); ++j) {
            std::istringstream ls(lines[j]);
            double idx, val;
            if (ls >> idx >> val)
                coefs.push_back(val);
        }

        orbitals.emplace_back(occ, coefs);
    }
    number_of_occupied_orbitals = n;

    std::cout << "Size of std::vector<Orbital> orbitals = " << orbitals.size() << std::endl;
}

double Molecule::scfp_density_at_point(const std::array<double, 3>& point) const {
    if (!density_matrix.empty()) {
        std::vector<double> values;
        values.reserve(basis_functions.size());
        for (const auto& bf : basis_functions) {
            auto value = bf.value_at_point(point);
            values.push_back(value);
        }

        // ρ = φᵀ P φ
        double res = 0.0;
        for (size_t i = 0; i < values.size(); ++i) {
            for (size_t j = 0; j < values.size(); ++j) {
                if (i < density_matrix.size() && j < density_matrix[i].size())
                    res += values[i] * density_matrix[i][j] * values[j];
            }
        }
        return res;
    } else if (!orbitals.empty()) {
        double res = 0.0;
        int i = 0;
        for (const auto& orb : orbitals){
            if (orb.get_occupancy() == 0) {
                break;
            }
            double val = orbital_value_at_point(point, i);
            res += orb.get_occupancy() * val * val;
            i++;
        }
        // for (size_t i = 0; i < orbitals.size(); ++i) {
        //     const auto& orb = orbitals[i];
        //     if (orb.get_occupancy() == 0) continue;
        //     double val = orbital_value_at_point(point, i);
        //     res += orb.get_occupancy() * val * val;
        // }
        return res;
    } else {
        return 0.0;
    }
}

double Molecule::orbital_value_at_point(const std::array<double, 3>& point, size_t number) const {
    if (number >= orbitals.size()) return 0.0;

    const auto& coeffs = orbitals[number].get_coefficients();

    if (basis_functions.size() != coeffs.size()) {
        throw std::runtime_error("Class Molecule. The number of basis functions is not equal to the number of coefficients");
    }

    double sum = 0.0;
    auto N = basis_functions.size();
    for (size_t i = 0; i < N; ++i) {
        // if (coeffs[i] == 0) continue;

        sum += basis_functions[i].value_at_point(point) * coeffs[i];
    }

    return sum;
}
