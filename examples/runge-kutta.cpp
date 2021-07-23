#include "runge-kutta.h"
#include <iostream>
#include <math.h>

using std::vector;
using std::cout;
using std::endl;

std::function<double(double, double, double)> rkm(
        std::function<double(double, double)> func,
        const vector<vector<double>> &matrix) {
    /*
       - func is the function on the right side of the ODE
       - matrix is the butcher tableau of a particular
           runge kutta method
       - rkm() returns a function which computes the difference
           between one stage and the next (y_{n + 1} - y_n)
    */

    return [func, matrix](double x, double y, double h) {
        // size of tableau, array of each k_i,
        // and variable representing Σ b_i * k_i
        int size = matrix.size();
        vector<double> arr;
        double sum {0};

        // for each k_i
        for (int m {0}; m < size - 1; m++) {
            double arg {y};
            for (int n {0}; n < m; n++) {
                // computing y arg for function call
                arg += h * arr[n] * matrix[m][n + 1];
            }

            // compute k_i, add to sum, and store
            double k = func(x + h * matrix[m][0], arg);
            sum += k * matrix[size - 1][m + 1];
            arr.push_back(k);
        }

        return h * sum;
    };
}

void example() {
    // example from the following link:
    // https://math.okstate.edu/people/yqwang/teaching/math4513_fall11/Notes/rungekutta.pdf

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
    auto func = [](double x, double y){return y - pow(x, 2) + 1;};
    auto diff = rkm(func, matrix);

    // initializing the variables necessary between stages
    double x {0},                        // initial x
           y {0.5},                      // initial y
           h {0.5},                      // step size
           stop {2},                     // estimated x
           num {round((stop - x) / h)};  // number of stages

    // loop through each stage, updating x and y
    for (int n = 0; n < num; n++) {
        y += diff(x, y, h);
        x += h;

        cout << "y(" << x << ") = "
             << y << endl << endl;
    }

    cout << "The answer is " << y
         << "." << endl;
}
