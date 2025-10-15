#ifndef ORBITAL_H
#define ORBITAL_H

#include <vector>
#include <string>

class Orbital {
public:
    Orbital(double occupancy, const std::vector<double>& coefficients);

    std::string repr() const;
    double get_occupancy() const { return occupancy_; } 
    std::vector<double> get_coefficients() const { return coefficients_; }

private:
    double occupancy_;
    std::vector<double> coefficients_;
};

#endif // ORBITAL_H
