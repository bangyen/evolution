#ifndef EVOLVE_H
#define EVOLVE_H

#include <vector>

const double PI = 3.14159265358979323846;

double alpha(double square, double lambda);
double step(int num, int total);
double value(int num, int total);
std::vector<double> evolve(int argc, char** argv);

#endif