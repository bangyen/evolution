#include "include/gpd.h"
#include "include/evolve.h"
#include <vector>
#include <string>
#include <functional>
#include <fstream>
//#include <sciplot/sciplot.hpp>

using std::vector;
using std::function;

extern vector<double> serial(
	double stop, double zeta, double t,
	function<double(double, double, double)> func,
	double c, double l, int num);

int main(int argc, char** argv) {
	int    num  = argc > 1 ? atoi(argv[1]) : 100;
	auto   func = gpdHuplus;
	double stop = 1;
	double zeta = 0.0001;
	double t    = -0.1;
	double c    = 4 / 3.0;
	double l    = 0.246;
	double a    = 0;
	double b    = 0;

	vector<double> val, init, res
		= serial(stop, zeta, t, func, c, l, num);

	std::fstream fs;
	fs.open ("output.csv", std::fstream::in | std::fstream::out | std::fstream::app);
	fs << "x\t" << "init\t" << "res\n";

	for (int k = 0; k < num; k++) {
		double temp = value(k, num);
		double h    = step(k - 1, num);

		val.push_back(temp);
		init.push_back(func(temp, zeta, t));

		a += init[k] * h;
		b +=  res[k] * h;

		fs << val[k] << "\t" << init[k] << "\t" << res[k] << "\n";
	}

	fs.close();
	/*
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
	*/
	return 0;
}
