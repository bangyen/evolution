#include "evolve.h"
#include <vector>
#include <cmath>

const double PI = 3.14159265358979323846;

double alpha(double square, double lambda) {
    return 4 * PI / (9 * log(square / pow(lambda, 2)));
}

// translation of alpha_sub.f
double alpha(double scale) {
    std::vector<double> xmas2, xmass = {1.4, 4.5, 170};
    int nflavor;

    for (auto k : xmass)
        xmas2.push_back(k * k);

    if (scale > xmas2[2])
        nflavor = 6;
    else if (scale > xmas2[1] && scale <= xmas2[2])
        nflavor = 5;
    else if (scale > xmas2[0] && scale <= xmas2[1])
        nflavor = 4;
    else
        nflavor = 3;

    double x, alphainv,
           xlam  = 0.215,
           beta0 = 11 - nflavor * 2 / 3.;

    switch (nflavor) {
        case 6: {
            double beta1 = 51 - nflavor * 19 / 3.;
            double x6 = log(xmass[2] / pow(xlam, 2));
            double add6 = 23 / 12. * x6 / (1 - 348 / 529. * log(x6) / x6)
                - 21 / 12. * x6 / (1 - 234 / 441. * log(x6) / x6);

            x = log(scale / pow(xlam, 2));
            alphainv = beta0 / 4 * x
                / (1 - 2 * beta1 / pow(beta0, 2) * log(x) / x) + add6;

            return 1 / alphainv;
        }
        case 5:
            xlam *= pow(xmass[1] / xlam, -2 / 23.);
            break;
        case 3:
            xlam *= pow(xmass[0] / xlam, 2 / 27.);
            break;
    }

    x = log(scale / pow(xlam, 2));
    alphainv = beta0 / 4 * x;

    return 1 / alphainv;
}

/*
   distance between current x value
     and the next
   
   num   = index of x value
   total = total number of x values
*/
double step(int num, int total) {
    int part = ++total / 3 * 2;

    if (++num < part)
        return 0.2 / (pow(2, part) - 1)
            * pow(2, num);
    return 0.8 / (total - part);
}

/*
   current x value

   num   = index of x value
   total = total number of x values
*/
double value(int num, int total) {
    int part = ++total / 3 * 2;

    if (++num < part)
        return 0.2 / (pow(2, part) - 1)
            * (pow(2, num) - 1);
    return 0.2 + 0.8 / (total - part)
        * (num - part);
}