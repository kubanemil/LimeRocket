from CEA_Wrap import Fuel, Oxidizer, RocketProblem, DataCollector
import matplotlib.pyplot as plt
from time import sleep
import os, math

keros = Fuel("RP-1", wt_percent=80)
ox = Oxidizer("O2(L)", temp=100, wt_percent=80)

chamber_pressure = 10
problem = RocketProblem(pressure=chamber_pressure, materials=[keros, ox], pressure_units="bar", o_f=1.0, pip=chamber_pressure)
results = problem.run()

exit_vel = results.son * results.mach
exit_pressure = results.p
exit_temp = results.t
chamber_pressure = results.c_p
chamber_temp = results.c_t
print(exit_vel, "m/s,", chamber_pressure,"Pa.", chamber_temp, "K.")
print(results)

def At_per_mdot(M, Tt, y, Pt):
    R = 8.31/M
    Pt = Pt*100000
    return math.sqrt((R*Tt)/(y/9.8)) / Pt

am = At_per_mdot(M=0.015, Tt=results.t_t, y=results.gamma, Pt=results.t_p)
r = math.sqrt((am*0.5)/3.14)
print(r*200)
# vs = []
# ofs = [of/10 for of in range(1, 30)]
# p = 10
# for of in ofs:
#     problem = RocketProblem(pressure=p, materials=[keros, ox], pressure_units="bar", o_f=of,
#                             pip=p)
#     results = problem.run()
#     vs.append(results.c_t)
#     for file in os.listdir():
#         if "my_output" in file:
#             os.remove(file)
#
# print(vs)
# plt.plot(ofs, vs)
# plt.show()
