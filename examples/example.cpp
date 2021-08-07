#include "../common/runge-kutta.h"
#include <iostream>
#include <chrono>
#include <math.h>

using std::vector;
using std::cout;
using std::endl;

int main() {
    // example from the following link:
    // https://math.okstate.edu/people/yqwang/teaching/math4513_fall11/Notes/rungekutta.pdf

    // the ODE function and resultant stage function
    auto func = [](double x, double y) {return y - pow(x, 2) + 1;};
    auto diff = rkm(func, rk4);

    // initializing the variables necessary between stages
    double x    {0},                      // initial x
           y    {0.5},                    // initial y
           h    {0.5},                    // step size
           stop {2},                      // estimated x
           num  {round((stop - x) / h)};  // number of stages

 // loop through each stage, updating x and y
    for (int n = 0; n < num; n++) {
        y += diff(x, y, h);
        x += h;

        cout << "y(" << x << ") = "
            << y << endl;
    }

    cout << "The answer is " << y
        << "." << endl;

    return 0;
}