#include "../common/gpd.h"
#include "evolve.h"
#include <vector>
#include <string>
#include <functional>
#include <sciplot/sciplot.hpp>

using std::vector;
using std::function;

extern vector<double> serial(
    vector<double> val,
    function<double(double, double, double)> func,
    double stop, double zeta, double t,
    double c, double l);

int main(int argc, char** argv) {
    int    num;
    auto   func = gpdHu;
    double a = 0, b = 0,
           stop = 1,
           zeta = 0,
           t    = 0,
           c    = 4 / 3.0,
           l    = 0.246;
    bool   xu   = true,
           csv  = true;
    
    vector<double> val, init, res;

    if (csv) {
        std::ifstream input("data.dat");
        std::istream& s = input;
        std::string str;

        while (std::getline(s, str)) {
            std::stringstream ss(str);
            if (ss.good())
            {
                std::string substr;
                std::getline(ss, substr, ' ');
                val.push_back(std::stod(substr));
            }
        }

        num = val.size();
    } else {
        num = argc > 1 ? std::stoi(argv[1]) : 100;
        for (int k = 0; k < num; k++)
            val.push_back(value(k, num));
    }

    res = serial(val, func, stop, zeta, t, c, l);
    val.push_back(1);

    for (int k = 0; k < num; k++) {
        double temp = val[k];
        double h    = val[k + 1] - temp;
        init.push_back(func(temp, zeta, t));

        if (xu)
            init[k] *= temp;
        else
            res[k] /= temp;

        a += init[k] * h;
        b +=  res[k] * h;
    }

    sciplot::Vec x(val.data(),  num),
                 y(init.data(), num),
                 z(res.data(),  num);
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
    plot.show();
    
    return 0;
}
