#include "BasisFunction.h"
#include <cmath>
#include <sstream>
#include <stdexcept>

#define R_THRESHOLD 10

BasisFunction::BasisFunction(const std::vector<double>& coefficients,
                             const std::vector<double>& exponents,
                             const std::string& shell,
                             const std::array<double, 3>& position,
                             const std::string& index,
                             const std::string& label)
    : coefficients_(coefficients),
      exponents_(exponents),
      shell_(shell),
      position_(position),
      index_(index),
      label_(label) {}

std::string BasisFunction::repr() const {
    std::ostringstream oss;
    oss << "Function label = " << label_ << "\n"
        << "position = (" << position_[0] << ", " << position_[1] << ", " << position_[2] << ")\n"
        << "coefficients = [";
    for (size_t i = 0; i < coefficients_.size(); ++i) {
        oss << coefficients_[i];
        if (i + 1 < coefficients_.size()) oss << ", ";
    }
    oss << "]\nexponents = [";
    for (size_t i = 0; i < exponents_.size(); ++i) {
        oss << exponents_[i];
        if (i + 1 < exponents_.size()) oss << ", ";
    }
    oss << "]\nshell = " << shell_ << "\nindex = " << index_;
    return oss.str();
}

double BasisFunction::value_at_point(const std::array<double, 3>& point) const {
    if (shell_ == "s") {
        return s(point);
    } else if (shell_ == "p") {
        return p(point);
    } else if (shell_ == "d") {
        return d(point);
    } else {
        throw std::runtime_error("Class BasisFunction: Unknown orbital shell");
    }
}

double BasisFunction::s(const std::array<double, 3>& point) const {
    double dx = position_[0] - point[0];
    double dy = position_[1] - point[1];
    double dz = position_[2] - point[2];
    double squared_radius_vector = dx * dx + dy * dy + dz * dz;

    // if (squared_radius_vector > R_THRESHOLD)
    //     return 0.0;

    // auto it = s_res_dict_.find(squared_radius_vector);
    // if (it != s_res_dict_.end())
    //     return it->second;

    double res = 0.0;
    for (size_t i = 0; i < coefficients_.size(); ++i)
        res += coefficients_[i] * std::exp(-exponents_[i] * squared_radius_vector);

    // s_res_dict_[squared_radius_vector] = res;
    return res;
}

double BasisFunction::p(const std::array<double, 3>& point) const {
    double dx = point[0] - position_[0];
    double dy = point[1] - position_[1];
    double dz = point[2] - position_[2];
    double squared_radius_vector = dx * dx + dy * dy + dz * dz;

    // if (squared_radius_vector > R_THRESHOLD)
    //     return 0.0;

    double res = 0.0;
    if (index_ == "x") {
        for (size_t i = 0; i < coefficients_.size(); ++i)
            res += coefficients_[i] * dx * std::exp(-exponents_[i] * squared_radius_vector);
    } else if (index_ == "y") {
        for (size_t i = 0; i < coefficients_.size(); ++i)
            res += coefficients_[i] * dy * std::exp(-exponents_[i] * squared_radius_vector);
    } else if (index_ == "z") {
        for (size_t i = 0; i < coefficients_.size(); ++i)
            res += coefficients_[i] * dz * std::exp(-exponents_[i] * squared_radius_vector);
    } else {
        throw std::runtime_error("Class BasisFunction: Unknown p-orbital index");
    }

    return res;
}

double BasisFunction::d(const std::array<double, 3>& point) const {
    // return 0.0;

    double dx = point[0] - position_[0];
    double dy = point[1] - position_[1];
    double dz = point[2] - position_[2];
    double squared_radius_vector = dx * dx + dy * dy + dz * dz;

    double res = 0.0;

    if (index_ == "z2") {
        for (size_t i = 0; i < coefficients_.size(); ++i)
            res += coefficients_[i] / std::sqrt(3.0) * (dz*dz - 0.5*dx*dx - 0.5*dy*dy) * std::exp(-exponents_[i] * squared_radius_vector);
    } else if (index_ == "xz") {
        for (size_t i = 0; i < coefficients_.size(); ++i)
            res += coefficients_[i] * dx * dz * std::exp(-exponents_[i] * squared_radius_vector);
    } else if (index_ == "yz") {
        for (size_t i = 0; i < coefficients_.size(); ++i)
            res += coefficients_[i] * dy * dz * std::exp(-exponents_[i] * squared_radius_vector);
    } else if (index_ == "x2y2") {
        for (size_t i = 0; i < coefficients_.size(); ++i)
            res += coefficients_[i] / 2 * (dx*dx - dy*dy) * std::exp(-exponents_[i] * squared_radius_vector);
    } else if (index_ == "xy") {
        for (size_t i = 0; i < coefficients_.size(); ++i)
            res += coefficients_[i] * dx * dy * std::exp(-exponents_[i] * squared_radius_vector);
    } else {
        throw std::runtime_error("Class BasisFunction: Unknown d-orbital index");
    }
    
    return res;
}
