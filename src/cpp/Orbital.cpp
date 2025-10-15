#include "Orbital.h"
#include <sstream>

Orbital::Orbital(double occupancy, const std::vector<double>& coefficients)
    : occupancy_(occupancy), coefficients_(coefficients) {}

std::string Orbital::repr() const {
    std::ostringstream oss;
    oss << "Occupancy = " << occupancy_ << "\ncoefficients = [";
    for (size_t i = 0; i < coefficients_.size(); ++i) {
        oss << coefficients_[i];
        if (i + 1 < coefficients_.size())
            oss << ", ";
    }
    oss << "]";
    return oss.str();
}
