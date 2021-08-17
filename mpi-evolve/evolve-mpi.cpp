#include "../common/runge-kutta.h"
#include "evolve.h"
#include <functional>
#include <vector>
#include <mpi.h>
#include <cmath>

using std::vector;
using std::function;

extern vector<double> evolve(
    double stop, double zeta, double t,
    function<double(double, double, double)> func,
    double c, double l) {
    // parallelization setup
    vector<double> res;
    int size, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    res.resize(size);

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
    double w = log(0.09362);
    double dw = (log(stop) - w) / size;
    double x = value(rank, size);
    double u = func(x, zeta, t);

    // the ODE function and resultant stage function
    auto stage = [&](double w2, double u2) {
        // each process collects u-values from each other
        MPI_Allgather(&u, 1, MPI_DOUBLE, res.data(), 1,
            MPI_DOUBLE, MPI_COMM_WORLD);

        double sum = u2 * (2 * log(1 - x) + 1.5);

        // compute riemann sum
        for (int k = rank; k < size; k++) {
            double y = value(k + 1, size);
            double temp = ((pow(x / y, 2) + 1)
                * res[k] - 2 * u2) / (y - x);
            sum += temp * step(k, size);
        }

        return c * alpha(exp(w2), l) * sum / (2 * PI);
    };
    auto diff = rkm(stage, rk4);

    // loop through each stage, updating w and u
    for (int n = 0; n < size; n++) {
        u += diff(w, u, dw);
        w += dw;
    }

    // root process collects each u and stores in res
    MPI_Gather(&u, 1, MPI_DOUBLE, res.data(), 1,
        MPI_DOUBLE, 0, MPI_COMM_WORLD);

    return res;
}