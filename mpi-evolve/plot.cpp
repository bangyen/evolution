#include "../common/gpd.h"
#include "evolve.h"
#include <vector>
#include <string>
#include <functional>
#include <sciplot/sciplot.hpp>

using std::vector;
using std::function;

extern vector<double> serial(
    double stop, double zeta, double t,
    function<double(double, double, double)> func,
    double c, double l, int num);

int main(int argc, char** argv) {
    std::ifstream file("gpd-u.csv");
    int    num  = argc > 1 ? std::stoi(argv[1]) : 100;
    auto   func = gpdHu;
    double a = 0, b = 0, d = 0, 
           stop = 1,
           zeta = 0.0001,
           t    = -0.1,
           c    = 4 / 3.0,
           l    = 0.246;
    bool   xu = false;

    vector<double> altx, alty, val, init, res
        = serial(stop, zeta, t, func, c, l, num);

    for (int k = 0; k < num; k++) {
        double temp = value(k, num);
        double h    = step(k - 1, num);

        val.push_back(temp);
        init.push_back(func(temp, zeta, t));

        if (xu) {
            init[k] *= temp;
            res[k] *= temp;
        }

        a += init[k] * h;
        b +=  res[k] * h;
    }

    while (file.good()) {
        std::string str;
        file >> str;
        
        int ind  = str.find(',');
        double w = std::stod(str.substr(0, ind)),
               u = std::stod(str.substr(ind + 1, str.length() - ind - 1));

        if (xu)
            u *= w;

        altx.push_back(w);
        alty.push_back(u);

        d += u * 0.01;
    }

    sciplot::Vec x(val.data(),  num),
                 y(init.data(), num),
                 z(res.data(),  num),
                 n(altx.data(), altx.size()),
                 m(alty.data(), alty.size());
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

    plot.drawCurve(x, y).label("Initial (Area: "   + std::to_string(a) + ")");
    plot.drawCurve(x, z).label("Estimated (Area: " + std::to_string(b) + ")");
    plot.drawCurve(n, m).label("CSV (Area: "       + std::to_string(d) + ")");
    plot.show();

    return 0;
}
