// c++ translation of fortran.f95
#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>

using std::string;
using std::cout;
using std::endl;
using std::setw;

// approximates dirac delta function
double delta(double x, double dx) {
    return x <= dx ? 1 / dx : 0;
}

// analytical solution
double phiex(double x) {
    return exp(-x) * sin(2 * x) / 2;
}

int main() {
    // set precision for output
    // and initialize data
    std::cout.precision(3);
    int    nstep   {2000};
    double L       {4};

    double dx      {L / nstep};
    int    kp      {(int) (nstep / 40.0)};
    double phi0    {0};
    double uprime0 {1};
    double phi1    {phi0 + dx * uprime0};

    // print column names and initial conditions
    cout << "  x          phi  "
         << "      phiex  ratio" << endl
         << "------------------"
         << "------------------" << endl
         << "  "     << 0
         << setw(13) << phi0
         << setw(13) << phiex(0)
         << "    "   << phi0 / phiex(0)
         << endl;

    for (int i = 2; i <= nstep; i++) {
        double x      {dx * i};
        double dphidx {(phi1 - phi0) / dx};
        double phi2;

        if (i == 1)
            phi2 = 2 * phi1 - phi0
                   + pow(dx, 2)
                   * (- 2 * dphidx
                      - 5 * phi1
                      + 1 / dx);
        else
            phi2 = 2 * phi1 - phi0
                   + pow(dx, 2)
                   * (- 2 * dphidx
                      - 5 * phi1);

        // for every kp steps, print data
        if (i % kp == 0) {
            cout.unsetf(std::ios::fixed);
            cout << std::noshowpos
                 << setw(3)  << x
                 << std::scientific
                 << setw(13) << phi2
                 << setw(13) << phiex(x)
                 << std::fixed
                 << "  "     << phi0 / phiex(x)
                 << endl;
        }

        phi0 = phi1;
        phi1 = phi2;
    }
}