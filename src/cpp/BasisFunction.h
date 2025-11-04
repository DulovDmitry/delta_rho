#ifndef BASISFUNCTION_H
#define BASISFUNCTION_H

#include <vector>
#include <string>
#include <unordered_map>
#include <array>

class BasisFunction {
public:
    BasisFunction(const std::vector<double>& coefficients,
                  const std::vector<double>& exponents,
                  const std::string& shell,
                  const std::array<double, 3>& position,
                  const std::string& index = "",
                  const std::string& label = "");
    ~BasisFunction() {};

    std::string repr() const;

    double value_at_point(const std::array<double, 3>& point) const;

private:
    std::vector<double> coefficients_;
    std::vector<double> exponents_;
    std::string shell_;
    std::array<double, 3> position_;
    std::string index_;
    std::string label_;

    std::unordered_map<double, double> s_res_dict_;

    double s(const std::array<double, 3>& point) const;
    double p(const std::array<double, 3>& point) const;
    double d(const std::array<double, 3>& point) const;
};

#endif // BASISFUNCTION_H
