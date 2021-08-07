#include <iostream>
#include <cmath>
#include <mpi.h>
#include "evolve.h"
#include "../common/gpd.h"
#include "../common/runge-kutta.h"

using std::vector;

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
    int part = total / 3 * 2;

    if (num < part)
        return 0.2 / (pow(2, part) - 1);
    return 0.8 / (total - part);
}

/*
   current x value

   num   = index of x value
   total = total number of x values
*/
double value(int num, int total) {
    int part = total / 3 * 2;

    if (num < part)
        return 0.2 / (pow(2, part) - 1)
            * (pow(2, num) - 1);
    return 0.2 + 0.8 / (total - part)
        * (num - part);
}

vector<double> evolve(int argc, char** argv) {
    // parallelization setup
    vector<double> res;
    int size, rank;

    MPI_Init(&argc, &argv);
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
    double c    = 4 / 3.0;
    double zeta = 0.0001;
    double t    = -0.1;
    double w    = log(0.09362);
    double stop = log(1);
    double l    = 0.246;
    double dw   = (stop - w) / size;
    double x    = value(rank, size);
    double u    = gpdHuplus(x, zeta, t);

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
        // each process collects u-values from each other
        MPI_Allgather(&u, 1, MPI_DOUBLE, res.data(), 1,
            MPI_DOUBLE, MPI_COMM_WORLD);

        double sum = u2 * (2 * log(1 - x) + 1.5);

        // compute riemann sum
        for (int k = rank; k < size; k++) {
            double y    = value(k + 1, size);
            double temp = ((pow(x / y, 2) + 1)
                * res[k] - 2 * u2) / (y - x);
            sum += temp * step(k, size);
        }

        return c * alpha(exp(w2), l) * sum / (2 * PI);
    };
    auto diff = rkm(func, matrix);

 // loop through each stage, updating w and u
    for (int n = 0; n < size; n++) {
        u += diff(w, u, dw);
        w += dw;
    }

    // root process collects each u and stores in res
    MPI_Gather(&u, 1, MPI_DOUBLE, res.data(), 1,
        MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0)
        for (int n = 1; n < size; n++)
            std::cout << "x: "   << value(n, size)
                      << ", u: " << res[n]
                      << std::endl;

    MPI_Finalize();
    return res;
}