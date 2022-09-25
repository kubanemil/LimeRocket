import matplotlib.pyplot as plt
import numpy as np
from math import e
from rockets import Rocket

width = 10
rocket = Rocket(tank_h=100, ox_is_liquid=True, lox_to_fuel=1, width=width)
dmdt = 0.0134 #kg/s
u = 800 #m/s
M_init = 0.0450 #kg
M_pt = 6.9333*0.01*M_init #kg
k = 0.00023

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
# plt.plot(dmdts, heights)
# plt.show()
