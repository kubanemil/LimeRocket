import matplotlib.pyplot as plt
import numpy as np
from math import e
from rockets import Rocket


width = 0.14 #meters
rocket = Rocket(width=width, others=2, fuel_height=0.9)
engine = rocket.Engine()
tanks = rocket.Tanks()
dmdt = engine.mdot
u = engine.exit_vel
M_init = rocket.total_mass()
M_pt = tanks.lox_mass()+tanks.fuel_mass() #kg
k = rocket.K()

print(M_init, M_pt, k)
def hc_vt_bt_th(dmdt=dmdt, u=u, k=k, M_init=M_init, M_pt=M_pt): #equation is valid
    Thrust = dmdt * u; btime = M_pt / dmdt
    Mb = M_init - (M_pt / 2); Mc = M_init - M_pt

    q = ((Thrust - (Mb * 9.8)) / k) ** (1 / 2)
    p = (2 * k * q) / Mb

    e_inx = e ** (-(p * btime))
    vt = q * (1 - e_inx) / (1 + e_inx)

    ln_inx = ((Mc * 9.8) + (k * (vt ** 2))) / (Mc * 9.8)
    hc = (Mc / (2 * k)) * np.log(ln_inx)
    return [hc, vt, btime, Thrust]

def g(height_km):
    return 9.8/(1+(2*height_km/6371))

dmdts = []
btimes = []
velocities = []
heights = []
us = []
# for t in range(M_pt+1, M_pt+30):
#     # print(hc_vt_bt_th(u=t))
#     dmdts.append(t)
#     heights.append(hc_vt_bt_th(M_init=t)[0])

print(hc_vt_bt_th())
print(hc_vt_bt_th(k=10**(-15)))
g = g(hc_vt_bt_th()[0]/1000)
print(g)