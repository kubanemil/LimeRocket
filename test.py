from math import sqrt, pi
from math import sin, cos, tan

pi = round(pi, 3)
def integral(func, rng):
    dx = 0.00001
    A = 0
    for x in range(int(rng[0]/dx), int(rng[1]/dx), 1):
        A += func(x*dx)*dx
    return A


def K(R, Cd=0.75):
    air_dens = 1.2  # kg/m3
    return (Cd*air_dens*pi/2)*(R**2)


def a_func(v):
    M = 10  # kg
    k = K(R=0.07, Cd=0.8)
    return -(k*(v**2)/M) - 9.8


def h_throw_up():
    dt = 0.001
    H = 0
    v = 10
    T = 0
    while True:
        a = a_func(v)
        v = v + a*dt
        H += v*dt
        T += dt
        if v < 0:
            return [H, T]

print(h_throw_up())

# second_order = lambda i: i**10
# g = integral(second_order, [0, 10])
# print(g)
# a_func=lambda t: t*10
# print(v(a_func, [0, 10]))
# print(cos(2*3.1415))
