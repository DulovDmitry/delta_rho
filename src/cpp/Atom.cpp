#include "Atom.h"
#include <sstream>
#include <stdexcept>

Atom::Atom(const std::array<double, 3>& position,
           const std::string& label,
           int number,
           double charge)
    : position_(position), label_(label), number_(number), charge_(charge) {}

std::string Atom::repr() const {
    std::ostringstream oss;
    oss << "Atom label = " << label_ << "\n"
        << "number = " << number_ << "\n"
        << "position = (" << position_[0] << ", "
        << position_[1] << ", " << position_[2] << ")";
    return oss.str();
}
