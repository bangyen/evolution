#include "../common/gpd.h"
#include "evolve.h"
#include <mpi.h>
#include <vector>
#include <string>
#include <functional>
#include <sciplot/sciplot.hpp>

using std::vector;
using std::function;

extern vector<double> evolve(
    double stop, double zeta, double t,
    function<double(double, double, double)> func,
    double c, double l);

int main(int argc, char** argv) {
    int size, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    auto   func = gpdHu;
    double a = 0, b = 0,
           stop = 1,
           zeta = 0,
           t    = 0,
           c    = 4 / 3.0,
           l    = 0.246;
    bool   xu   = false;

    vector<double> val, init, res
        = evolve(stop, zeta, t, func, c, l);

    if (rank == 0) {
        for (int k = 0; k < size; k++) {
            double temp = value(k, size);
            double h = step(k - 1, size);

            val.push_back(temp);
            init.push_back(func(temp, zeta, t));

            if (xu) {
                init[k] *= temp;
                res[k] *= temp;
            }

            a += init[k] * h;
            b += res[k] * h;
        }

        sciplot::Vec x(val.data(), size),
            y(init.data(), size),
            z(res.data(), size);
        sciplot::Plot plot;

        plot.size(1440, 800);
        plot.xlabel("x");
        plot.ylabel("u(x, Q^2)");
        plot.xrange(0.0, 1.0);
        plot.yrange(0.0, 4.0);

        plot.legend()
            .atOutsideBottom()
            .fontSize(10)
            .displayHorizontal()
            .displayExpandWidthBy(2);

        plot.drawCurve(x, y).label("Initial (Area: " + std::to_string(a) + ")");
        plot.drawCurve(x, z).label("Estimated (Area: " + std::to_string(b) + ")");
        plot.save("evolve.pdf");
    }

    MPI_Finalize();
    return 0;
}
