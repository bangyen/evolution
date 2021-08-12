#ifndef RUNGE_KUTTA_H
#define RUNGE_KUTTA_H

#include <functional>
#include <vector>

/*
   the butcher tableau of RK4,
     leftmost column is nodes (c_i)
     and bottom row is weights (b_i)
     whereas the rest is the runge kutta matrix
*/
const std::vector<std::vector<double>> rk4 = {
    {0},
    {0.5, 0.5},
    {0.5, 0, 0.5},
    {1, 0, 0, 1},
    {0, 1. / 6, 1. / 3, 1. / 3, 1. / 6}
};

std::function<double(double, double, double)> rkm(
    std::function<double(double, double)> func,
    const std::vector<std::vector<double>> &matrix);

#endif