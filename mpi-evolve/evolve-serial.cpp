#include <iostream>
#include <cmath>
#include "evolve.h"
#include "../common/gpd.h"
#include "../common/runge-kutta.h"

using std::vector;

vector<double> serial(int num) {
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
    double c    = 4 / 3.0;
    double zeta = 0.0001;
    double t    = -0.1;
    double w    = log(0.09362);
    double stop = log(1);
    double l    = 0.246;
    double dw   = (stop - w) / num;

    /*
       initialize res with initial
       u-values for each x value
    */
    for (int n = 0; n < num; n++) {
        double x = value(n, num);
        temp.push_back(gpdHuplus(x, zeta, t));
    }

    /*
       the butcher tableau of RK4,
         leftmost column is nodes (c_i)
         and bottom row is weights (b_i)
         whereas the rest is the runge kutta matrix
    */
    vector<vector<double>> matrix = {
        {0},
        {0.5, 0.5},
        {0.5, 0, 0.5},
        {1, 0, 0, 1},
        {0, 1 / 6.0, 1 / 3.0, 1 / 3.0, 1 / 6.0}
    };

    // the ODE function and resultant stage function
    auto func = [&](double w2, double u2) {
        double x   = value(p, num);
        double sum = u2 * (2 * log(1 - x) + 1.5);

        // compute riemann sum
        for (int k = p; k < num; k++) {
            double y   = value(k + 1, num);
            double val = ((pow(x / y, 2) + 1)
                * res[k] - 2 * u2) / (y - x);
            sum += val * step(k, num);
        }

        return c * alpha(exp(w2), l) * sum / (2 * PI);
    };
    auto diff = rkm(func, matrix);

    // loop through each stage, updating w and u
    for (int n = 0; n < num; n++) {
        res = temp;

        for (p = 0; p < num; p++)
            temp[p] += diff(w, temp[p], dw);

        w += dw;
    }

    return temp;
}