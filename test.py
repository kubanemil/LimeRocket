from math import sqrt
import matplotlib.pyplot as plt

def area(r):
    return 3.14*(r**2)

def area_to_radius(A):
    return sqrt(A/3.14)

def mdot(Cd=0.8, r=0.03, Pt=280000, R=8.31, M=0.032, Tt=1280, y=1.2):
    ys = (2/(y+1))**((y+1)/(y-1))
    sqrt_ix = (2*M*ys)/(R*Tt)
    return Cd*area(r)*Pt*sqrt(sqrt_ix)

Ar = area(0.06)
print(Ar)
print(1/Ar)
r1 = [area_to_radius(16.589), area_to_radius(1), area_to_radius(3.8141)]
ix = 3/r1[0]
rs = [r*ix for r in r1]
mrs = [-r for r in rs]
print(rs)

l=5
xs = [0, l, l+2]
plt.plot(xs, rs, color='k')
plt.plot(xs, mrs, color='k')
plt.legend(loc="best")
plt.show()

