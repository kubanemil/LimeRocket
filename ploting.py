import matplotlib.pyplot as plt
import numpy as np
"""
u = speed of exhaust
Mt = mass of fuel/LOX
Mr = mass of rocket without fuel/LOX
t = time of engines work
"""
dmdt = 0.2 #kg/s
u = 2500 #m/s
M_init = 23 #kg
M_pt = 12 #kg
Thrust = dmdt*u
btime = M_pt/dmdt
print(dmdt, u, M_init, M_pt, Thrust, btime)

def v(t=btime):
    mass_inx = M_init/(M_init - (dmdt*t))
    return (  u*( np.log(mass_inx) )  ) - (9.8*t)

def H(v, t=btime):
    ho = v*t/2
    h1 = (v*v)/(2*9.8)
    return ho+h1


print("Bornout velocity:", v(), "m/s")
print("Coast height:", H(v()), "meters")

mrs = np.array([i for i in range(5, 100)])

vs = np.array([mr*mr for mr in mrs])
# Hs = np.array([H(V, t) for V in vs])
plt.plot(mrs, vs)
plt.show()
# print(Hs[-1])
# for i in range(1, 20):
#     v1 = v(u, Mt, Mr, t=i)
#     h1 = H(v1, t)
#     print(v1//1, h1//1)


