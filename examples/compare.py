# this file compares the accuracy
# and speed of using RK4 vs the code
# located at http://www.uprh.edu/rbaretti/IntegroDE13oct2014.htm
# for the example given

import matplotlib.pyplot as plot
from timeit import timeit
import numpy as np
import math
import sys


# analytic solution
def analytic(size, rep):
    result = []

    for k in range(rep):
        x = k * size
        result.append(
            math.exp(-x) / 2
            * math.sin(2 * x))

    return result


# python translation of difference.f95
def difference(size, rep):
    result = [0, size]
    phi_0 = 0
    phi_1 = size

    for _ in range(2, rep):
        prime = (phi_1 - phi_0) / size
        phi_2 = (2 * phi_1 - phi_0 + size ** 2
                 * (-2 * prime - 5 * phi_1))

        result.append(phi_2)
        phi_0 = phi_1
        phi_1 = phi_2

    return result


#  Equation: u' + 2u + 5∫_0^x u = H(x), u(0) = 0  #
#  Rewritten as second-order ODE:                 #
#    q'' + 2q' + 5q = H(x), q(0) = 0, q'(0) = 0   #
#  Further rewritten as two first-order ODEs:     #
#    y_1 = q  --→  y_1' = q' = y_2                #
#    y_2 = q'      y_2' = H(x) - 5q - 2q'         #
#                       = H(x) - 5y_1 - 2y_2      #
#    y_1(0) = 0, y_2(0) = 0                       #

def runge(h, rep):
    # b is a 1d array of size 2
    # representing [y_1, y_2]
    def f(a, b):
        res = -5 * b[0] - 2 * b[1]
        if a >= 0:
            res += 1
        return np.array([b[1], res])

    x = 0
    y = np.array([0.0, 0.0])
    result = [0]

    # runge() uses RK4 but with vectors
    for _ in range(1, rep):
        k_1 = h * f(x, y)
        k_2 = h * f(x + h / 2, y + k_1 / 2)
        k_3 = h * f(x + h / 2, y + k_2 / 2)
        k_4 = h * f(x + h, y + k_3)

        y += (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6
        result.append(y[1])
        x += h

    return result


if __name__ == '__main__':
    step = 0.002
    num = 3000

    if len(sys.argv) == 3:
        try:
            step = float(sys.argv[1])
            num = int(sys.argv[2])
        except ValueError:
            pass

    # get data from each function as
    # well as corresponding x-values
    x_list = [k * step for k in range(num)]
    a_arr = np.array(analytic(step, num))
    d_arr = np.array(difference(step, num))
    r_arr = np.array(runge(step, num))

    # determine the absolute
    # error of each function
    d_arr = [*map(abs, a_arr - d_arr)]
    r_arr = [*map(abs, a_arr - r_arr)]

    figure, axis_1 = plot.subplots()
    color = 'tab:red'
    axis_1.set_ylabel('difference', color=color,
                      labelpad=5)
    axis_1.tick_params(axis='y', labelcolor=color)
    axis_1.plot(x_list, d_arr, color=color)
    axis_1.set_title('Equation Solution Error (Absolute)',
                     fontsize=16, pad=15)

    axis_2 = axis_1.twinx()
    color = 'tab:blue'
    axis_2.set_ylabel('runge-kutta', color=color,
                      rotation=270, labelpad=15)
    axis_2.tick_params(axis='y', labelcolor=color)
    axis_2.plot(x_list, r_arr, color=color)

    figure.tight_layout()
    plot.show()

    print('Difference: ', timeit(
        lambda: difference(step, num), number=5))
    print('Runge-Kutta:', timeit(
        lambda: runge(step, num), number=5))
