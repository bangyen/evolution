#include "runge-kutta.h"
#include <iostream>
#include <chrono>
#include <math.h>

using std::vector;
using std::cout;
using std::endl;

int main() {
    using sec   = std::chrono::duration<double>;
    using clock = std::chrono::system_clock;
    auto before = clock::now();
    vector<double> res;

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

    for (int k = 0; k < 10; k++) {
        // the ODE function and resultant stage function
        auto func = [k](double x, double y){return k * y;};
        auto diff = rkm(func, matrix);

        // initializing the variables necessary between stages
        double x {0},                        // initial x
               y {1},                        // initial y
               h {0.000005},                 // step size
               stop {2},                     // estimated x
               num {round((stop - x) / h)};  // number of stages

        // loop through each stage, updating x and y
        for (int n = 0; n < num; n++) {
            y += diff(x, y, h);
            x += h;
        }

        res.push_back(y);
    }

    for (int k = 0; k < 10; k++)
        cout << "e ^ " << 2 * k << " = "
             << res[k] << endl;

    sec time = clock::now() - before;
    cout << "Duration: "
         << time.count()
         << "s" << endl;

    return 0;
}