#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <array>

class Atom {
public:
    Atom(const std::array<double, 3>& position,
         const std::string& label = "",
         int number = 0,
         double charge = 0.0);

    std::string repr() const;

    // Геттеры
    const std::array<double, 3>& get_position() const { return position_; }
    const std::string& get_label() const { return label_; }
    int get_number() const { return number_; }
    double get_charge() const { return charge_; }

private:
    std::array<double, 3> position_;
    std::string label_;
    int number_;
    double charge_;
};

#endif // ATOM_H
