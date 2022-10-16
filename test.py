from math import sqrt, pi
from math import sin, cos, tan
from rockets import Rocket


pi = round(pi, 3)
def integral(func, rng):
    dx = 0.00001
    A = 0
    for x in range(int(rng[0]/dx), int(rng[1]/dx), 1):
        A += func(x*dx)*dx
    return A




print(flight())

# second_order = lambda i: i**10
# g = integral(second_order, [0, 10])
# print(g)
# a_func=lambda t: t*10
# print(v(a_func, [0, 10]))
# print(cos(2*3.1415))
