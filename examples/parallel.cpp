#include <mpi.h>
#include <math.h>
#include <chrono>
#include <iostream>
#include "./runge-kutta.h"

using std::vector;
using std::cout;
using std::endl;

int main(int argc, char **argv) {
    // rename clock and get starting time
    using sec   = std::chrono::duration<double>;
    using clock = std::chrono::system_clock;
    auto before = clock::now();

    // parallelization setup
    vector<double> res;
    int size, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
        res.resize(size);

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
        {0, 1/6.0, 1/3.0, 1/3.0, 1/6.0}
    };

    // the ODE function and resultant stage function
    auto func = [rank](double x, double y){return rank * y;};
    auto diff = rkm(func, matrix);

    // initializing the variables necessary between stages
    double x    {0},                      // initial x
           y    {1},                      // initial y
           h    {0.000005},               // step size
           stop {2},                      // estimated x
           num  {round((stop - x) / h)};  // number of stages

    // loop through each stage, updating x and y
    for (int n = 0; n < num; n++) {
        y += diff(x, y, h);
        x += h;
    }

    // root process collects each y and stores in res
    MPI_Gather(&y, 1, MPI_DOUBLE, res.data(), 1,
               MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // root process prints out result and time
    if (rank == 0) {
        for (int k = 0; k < size; k++)
            cout << "e ^ " << 2 * k << " = "
                 << res[k] << endl;

        sec time = clock::now() - before;
        cout << "Duration: "
             << time.count()
             << "s" << endl;
    }

    MPI_Finalize();
    return 0;
}
