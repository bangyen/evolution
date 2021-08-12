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
	int    num  = argc > 1 ? atoi(argv[1]) : 100;
	auto   func = gpdHu;
	double stop = 1;
	double zeta = 0.0001;
	double t    = -0.1;
	double c    = 4 / 3.0;
	double l    = 0.246;
	double a = 0, b = 0;

	vector<double> val, init, res
		= serial(stop, zeta, t, func, c, l, num);

	for (int k = 0; k < num; k++) {
		double temp = value(k, num);
		double h    = step(k - 1, num);

		val.push_back(temp);
		init.push_back(func(temp, zeta, t) * temp);
		res[k] = res[k] * temp;

		a += init[k] * h;
		b +=  res[k] * h;
	}

	sciplot::Vec x(val.data(),  num),
				 y(init.data(), num),
				 z(res.data(),  num);
	sciplot::Plot plot;

	plot.size(720, 400);
	plot.xlabel("x");
	plot.ylabel("u(x, Q^2)");
	plot.xrange(0.0, 1.0);
	plot.yrange(0.0, 4.0);

	plot.legend()
		.atTop()
		.fontSize(10)
		.displayHorizontal()
		.displayExpandWidthBy(2);

	plot.drawCurve(x, y).label("Initial (Area: "   + std::to_string(a) + ")");
	plot.drawCurve(x, z).label("Estimated (Area: " + std::to_string(b) + ")");
	plot.show();

	return 0;
}
