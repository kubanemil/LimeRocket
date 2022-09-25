import matplotlib.pyplot as plt
import numpy as np

#kerosene: dynamic viscosity (w) = 0.00164 Pa*s     d=800 kg/m3
#water: w = 0.0010518 Pas*s         d=1000 kg/m3
R = 0.0025 #meters
L = 0.01 #meters
p0 = 6 * (10**5) #Pa
pf = 5 * (10**5) #Pa
w = 0.0010518 #Pa*s
d = 1000 #kg/m3 density of kerosene
def mdot(R=R, L=L, p0=p0, pf=pf): #mass rate flow
    return (3.14*d*(R**4)*(p0-pf))/(8*w*L)

def p_diff(mdot, L=L, R=R):
    return (8*w*L*mdot)/(3.14*d*(R**4))

# def P_diff(mdot, Cd, R, dens):
#     A = 3.14*(R**2)
#     sqrt_ix = (mdot/(Cd*A))**2
#     return sqrt_ix/(2*9.8*dens)


def v(Cd, p_diff, density):
    return Cd* np.sqrt(2*g*p_diff/density)

print(p_diff(mdot=0.1, L=1.5, R=0.0015))
print(p_diff(mdot=0.03,  L=1, R=0.00125))
print(mdot())

