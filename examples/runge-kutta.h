#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H

#include <functional>
#include <vector>

std::function<double(double, double, double)> rkm(
    std::function<double(double, double)> func,
    const std::vector<std::vector<double>> &matrix);
void example();

#endif