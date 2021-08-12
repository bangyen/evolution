#include "../include/evolve.h"
#include <cmath>

const double PI = 3.14159265358979323846;

double alpha(double square, double lambda) {
    return 4 * PI / (9 * log(square / pow(lambda, 2)));
}

/*
   distance between current x value
     and the next
   
   num   = index of x value
   total = total number of x values
*/
double step(int num, int total) {
    int part = ++total / 3 * 2;

    if (++num < part)
        return 0.2 / (pow(2, part) - 1)
            * pow(2, num);
    return 0.8 / (total - part);
}

/*
   current x value

   num   = index of x value
   total = total number of x values
*/
double value(int num, int total) {
    int part = ++total / 3 * 2;

    if (++num < part)
        return 0.2 / (pow(2, part) - 1)
            * (pow(2, num) - 1);
    return 0.2 + 0.8 / (total - part)
        * (num - part);
}
