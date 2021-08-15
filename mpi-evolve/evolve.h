#ifndef EVOLVE_H
#define EVOLVE_H

extern const double PI;

double alpha(double square, double lambda);
double alpha(double scale);
double step(int num, int total);
double value(int num, int total);

#endif