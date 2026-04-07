#include <iostream>
#include src.hpp

const int interface_size = 3;
const int connection_size = 3;
int fromA[connection_size] = {1, 1, 2};
int toA[connection_size] = {2, 3, 3};
fraction resistanceA[connection_size] = {fraction(1, 2), fraction(1, 4), fraction(2)};
fraction currentA[interface_size] = {fraction(2), fraction(1), fraction(-3)};
fraction voltageA[interface_size] = {fraction(1), fraction(2), fraction(1, 2)};

int main() {
    resistive_network instance = resistive_network(interface_size, connection_size, fromA, toA, resistanceA);
    std::cout << instance.get_equivalent_resistance(1, 2) << n; // 9/22
    std::cout << instance.get_voltage(2, currentA) << n;       // 10/11
    std::cout << instance.get_power(voltageA) << n;             // 33/8
}
