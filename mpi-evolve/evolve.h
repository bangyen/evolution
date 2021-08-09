#ifndef EVOLVE_H
#define EVOLVE_H

#include <functional>
#include <vector>

const double PI = 3.14159265358979323846;

double alpha(double square, double lambda);
double step(int num, int total);
double value(int num, int total);
std::vector<double> evolve(
    double stop, double zeta, double t,
    std::function<double(double, double, double)> func,
    double c, double l, int argc, char** argv);

#endif