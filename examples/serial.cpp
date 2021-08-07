#include "../common/runge-kutta.h"
#include <iostream>
#include <chrono>
#include <math.h>

using std::vector;
using std::cout;
using std::endl;

int main(int argc, char **argv) {
    /*
       rename clock, get starting time,
         set number of iterations
    */
    int val     = argc > 1 ? atoi(argv[1]) : 10;
    using sec   = std::chrono::duration<double>;
    using clock = std::chrono::system_clock;
    auto before = clock::now();
    vector<double> res;

    for (int k = 0; k < val; k++) {
        // the ODE function and resultant stage function
        auto func = [k](double x, double y){return k * y;};
        auto diff = rkm(func, rk4);

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

    // print values
    for (int k = 0; k < val; k++)
        cout << "e ^ " << 2 * k << " = "
             << res[k] << endl;

    // print time
    sec time = clock::now() - before;
    cout << "Duration: "
         << time.count()
         << "s" << endl;

    return 0;
}