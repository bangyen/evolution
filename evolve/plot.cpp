#include "include/gpd.h"
#include "include/evolve.h"
#include <vector>
#include <string>
#include <functional>
#include <fstream>

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
    double a    = 0,
		   b    = 0,
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
    }
    else {
        num = argc > 1 ? std::stoi(argv[1]) : 100;
        for (int k = 0; k < num; k++)
            val.push_back(value(k, num));
    }

    res = serial(val, func, stop, zeta, t, c, l);
    val.push_back(1);

	std::fstream fs;
	fs.open ("output.csv", std::fstream::in | std::fstream::out | std::fstream::app);
	fs << "x\t" << "init\t" << "res\n";

	for (int k = 0; k < num; k++) {
		double temp = val[k];
		double h = val[k + 1] - temp;
		init.push_back(func(temp, zeta, t));

		if (xu)
			init[k] *= temp;
		else
			res[k] /= temp;

		a += init[k] * h;
		b += res[k] * h;

		fs << val[k] << "\t" << init[k] << "\t" << res[k] << "\n";
	}

	fs.close();
	return 0;
}
