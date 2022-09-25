import matplotlib.pyplot as plt
import numpy as np

V = 0.001 #m3
M = 0.032 #kg/mole
R = 8.31 #J/mole.K
T = 300 #K
X = (R*T)/(M*V)
Po = 2000000
Mo = Po/X
print(Mo)
dmdt = 0.1
def P(t):
    return (Mo - (dmdt*t)) * X


ts = np.array([t/10 for t in range(0, 33)])
Ps = np.array([P(t)/100000 for t in ts])
plt.plot(ts, Ps)
plt.show()


def dPdt(Po, Vo, mdot, dens, t):
    nRTmdot = Po*Vo*mdot
    V_ix = (Vo + (mdot*t/dens))**2
    return - (nRTmdot/(dens*V_ix))