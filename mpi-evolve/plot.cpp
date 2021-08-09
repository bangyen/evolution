#include <sciplot/sciplot.hpp>
#include "../common/gpd.h"
#include "evolve.h"
#include <vector>
#include <functional>

using std::vector;
using std::function;

vector<double> serial(
	double stop, double zeta, double t,
	function<double(double, double, double)> func,
	double c, double l, int num);

int main(int argc, char** argv) {
	int    num  = argc > 1 ? atoi(argv[1]) : 50;
	auto   func = gpdHuplus;
	double stop = 1;
	double zeta = 0.0001;
	double t    = -0.1;
	double c    = 4 / 3.0;
	double l    = 0.246;

	vector<double> val, init, res
		= serial(stop, zeta, t, func, c, l, num);
	res = vector<double>(res.begin() + 1, res.end());

	for (int k = 1; k < num; k++) {
		double temp = value(k, num);

		val.push_back(temp);
		init.push_back(func(temp, zeta, t));
	}

	sciplot::Vec x(val.data(),  num - 1),
				 y(init.data(), num - 1),
				 z(res.data(),  num - 1);
	sciplot::Plot plot;

	plot.size(720, 400);
	plot.xlabel("x");
	plot.ylabel("u(x, Q^2)");
	plot.xrange(0.0, 1.0);
	plot.yrange(0.0, 5.0);

	plot.legend()
		.atTop()
		.fontSize(10)
		.displayHorizontal()
		.displayExpandWidthBy(2);

	plot.drawCurve(x, y).label("initial");
	plot.drawCurve(x, z).label("estimated");
	plot.show();

	return 0;
}
