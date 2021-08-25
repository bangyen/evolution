#include "../common/runge-kutta.h"
#include "evolve.h"
#include <iostream>
#include <cmath>

using std::vector;
using std::function;

extern vector<double> serial(
        vector<double> vec,
        function<double(double, double, double)> func,
        double stop, double zeta, double t,
        double c, double l) {
    vector<double> res, temp;
    int p;

    /*
       c    = color factor
       zeta = xbj
       w    = ln(Q ^ 2)
       stop = estimated w
       l    = lambda (for alpha eq)
       dw   = step size for w
       x    = x value for current process
       u    = estimated u(x, Q ^ 2)
    */
    int    num  = vec.size();
    double w    = log(0.09362);
    double dw   = (log(stop) - w) / num;

    /*
       initialize res with initial
       u-values for each x value
    */
    for (auto k : vec)
        temp.push_back(k * func(k, zeta, t));
    vec.push_back(1);

    // the ODE function and resultant stage function
    auto stage = [&](double w2, double u2) {
        double x   = vec[p];
        double sum = u2 * (2 * log(1 - x) + 1.5);

        // compute riemann sum
        for (int k = p; k < num; k++) {
            double y   = vec[k + 1];
            double val = ((pow(x / y, 2) + 1)
                * res[k] - 2 * u2) / (y - x);
            sum += val * (y - vec[k]);
        }

        return c * alpha(exp(w2)) * sum / (2 * PI);
    };
    auto diff = rkm(stage, rk4);

    // loop through each stage, updating w and u
    for (int n = 0; n < num; n++) {
        res = temp;

        for (p = 0; p < num; p++)
            temp[p] += diff(w, temp[p], dw);

        w += dw;
    }

    return temp;
}