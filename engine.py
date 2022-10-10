"""
The class model of rocket engine that returns CEA results, also
estimating the sizing of the engine.
"""
# TODO: Create a method to draw a rocket scheme with valid proportions
# TODO: Estimate how fast will the rocket reach to it's critical point of heating.

from tools import *
from CEA_Wrap import Fuel, Oxidizer, RocketProblem
import matplotlib.pyplot as plt
import numpy as np


#CEA verified. #Sizing verified.
class Engine:
    def __init__(self, ch_pressure=5, mdot=0.5, ch_radius=0.06, o_f=1):
        self.ch_pressure = ch_pressure  # bar
        self.mdot = mdot    # kg/s
        self.o_f = o_f    # oxidizer to fuel ratio
        self.ch_radius = ch_radius  # m
        self.ch_area = R_to_A(self.ch_radius)   # m2
        self.ch_length = 0.1
        self.thickness = 4 * (10 ** (-3))  # meters
        self.a_exit = 12   # degree
        self.a_throat = 30    # degree
        self.insul_thick = 0.2  # cm
        self.create_properties()
        self.create_engine_sizes()

    def results(self):
        results = RocketProblem(pressure=self.ch_pressure, materials=[Fuel('RP-1'), Oxidizer('O2(L)', temp=100)], \
                                pressure_units="bar", o_f=self.o_f, pip=self.ch_pressure).run()
        results.set_fac_ma = self.mdot / self.ch_area  # kg/s / m2
        for el in list(results.keys())[:12]:
            del results[el]
        return results


    def total_mass(self):
        total_surarea = self.ch_surarea + self.th_nozzle_surarea + self.exit_nozzle_surarea
        dens = 7900     # kg/m3
        t = self.thickness  # meters
        return round(total_surarea*t*dens, 3)

    def plot_engine(self, ix=999, save=False):
        lengths_cm = self.lengths*100; radiuses_cm = self.radiuses*100
        plt.rcParams["figure.figsize"] = [13.9, 6.5]
        plt.rcParams["figure.autolayout"] = True
        plt.figure(ix)
        plt.axis([0, 30, -8, 8])
        plt.axis('scaled')
        plt.plot(lengths_cm, radiuses_cm, "red")
        plt.plot(lengths_cm, -radiuses_cm, "red")
        plt.plot(lengths_cm, radiuses_cm+self.insul_thick, "black")
        plt.plot(lengths_cm, -(radiuses_cm+self.insul_thick), "black")
        plt.figtext(0.13, 0.93, s=f"Thrust={round(self.thrust/9.8, 2)} kg", fontsize='12')
        plt.figtext(0.13, 0.9, s=f"Engine's mass={round(self.total_mass(), 2)} kg", fontsize='12')
        plt.figtext(0.30, 0.93, s=f"Mass flow={round(self.mdot, 2)} kg/s", fontsize='12')
        plt.figtext(0.30, 0.9, s=f"O/F ratio={round(self.o_f, 2)}", fontsize='12')
        plt.figtext(0.45, 0.93, s=f"ISP={round(self.isp, 2)} sec",  fontsize='12')
        plt.figtext(0.45, 0.9, s=f"Exit velocity={round(self.exit_vel, 2)} m/s", fontsize='12')

        plt.figtext(0.24, 0.53, s=f"{round(self.ch_pressure, 2)} bar", fontsize='11')
        plt.figtext(0.24, 0.5, s=f"{round(self.ch_temp, 2)} K", fontsize='11')
        plt.figtext(lengths_cm[-3]/32, 0.53, s=f"{round(self.th_pressure, 2)} bar", fontsize='11')
        plt.figtext(lengths_cm[-3]/32, 0.5, s=f"{round(self.th_temp, 2)} K", fontsize='11')
        plt.figtext(lengths_cm[-1]/32, 0.53, s=f"{round(self.exit_pressure, 2)} bar", fontsize='11')
        plt.figtext(lengths_cm[-1]/32, 0.5, s=f"{round(self.exit_temp, 2)} K", fontsize='11')
        plt.xlabel("Length, cm")
        plt.ylabel("Radius, cm")
        if save:
            plt.savefig(f"engine_plots\\engine{ix}.png")
        else:
            plt.show()

    def __repr__(self):
        return f"Thrust: {round(self.thrust/9.8, 2)} kg.\n" + \
               f"Exit velocity: {round(self.exit_vel, 2)} m/s.\n" + \
            f"Chamber pressure and temperature: {self.ch_pressure} bar, {self.ch_temp} K.\n" + \
            f"Specific impulse, sea-level: {round(self.isp, 2)} seconds."

    def create_properties(self):
        results = self.results()
        self.exit_vel = results.son * results.mach;
        self.exit_son = results.son;
        self.exit_mach = results.mach;
        self.exit_pressure = results.p;
        self.exit_temp = results.t;
        self.exit_dens = results.rho
        self.ch_pressure = results.c_p;
        self.ch_temp = results.c_t;
        self.ch_dens = results.c_rho
        self.th_pressure = results.t_p;
        self.th_temp = results.t_t;
        self.th_dens = results.t_rho
        self.gamma = results.gamma;
        self.molar_mass = results.mw / 1000;
        self.isp = results.isp
        self.thrust = self.exit_vel * self.mdot

    def nozzle_surface_area(self, th_radius, R, a):
        tana = tan(a * pi / 180)  # degree divergence angle
        return (pi / tana) * ((R ** 2) - (th_radius ** 2))

    def create_engine_sizes(self):
        self.Ae_At = self.Ae_At()
        self.th_area = self.th_area()
        self.exit_area = self.Ae_At * self.th_area
        self.exit_radius = A_to_R(self.exit_area)
        self.th_radius = A_to_R(self.th_area)


        self.exit_length = (self.exit_radius - self.th_radius) / (tan(self.a_exit * pi / 180))
        self.th_length = (self.ch_radius - self.th_radius) / (tan(self.a_throat * pi / 180))

        self.exit_nozzle_surarea = self.nozzle_surface_area(self.th_radius, self.exit_radius, a=self.a_exit)
        self.th_nozzle_surarea = self.nozzle_surface_area(self.th_radius, self.ch_radius, a=self.a_throat)
        self.ch_surarea = (2 * pi * self.ch_radius * self.ch_length) + (pi * (self.ch_radius ** 2))

        self.radiuses = np.array([self.ch_radius, self.ch_radius, self.th_radius, self.exit_radius, 0])
        self.lengths = np.array(scaled_list([0, self.ch_length, self.th_length, self.exit_length, 0]))
        self.reversed_radiuses = np.array([self.exit_radius, self.th_radius, self.ch_radius, self.ch_radius, 0])
        self.reversed_lengths = np.array(scaled_list([0, self.exit_length, self.th_length, self.ch_length, 0]))
        self.total_length = self.ch_length + self.th_length + self.exit_length

    def Ae_At(self):
        k = self.gamma
        pip = 1 / self.ch_pressure
        first_ix = ((k + 1) / 2) ** (1 / (k - 1)) * pip ** (1 / k)
        sqrt_ix = (k + 1) / (k - 1) * (1 - pip ** ((k - 1) / k))
        At_Ae = first_ix * (sqrt_ix) ** 0.5
        return 1 / At_Ae

    def th_area(self):
        Pc = self.ch_pressure * 100000;
        y = self.gamma
        Tc = self.ch_temp;
        mdot = self.mdot;
        R = 8.314 / self.molar_mass
        first = mdot / Pc;
        second = (R * Tc) / y;
        third = (2 / (y + 1));
        fourth = (y + 1) / (y - 1)
        return first * (second / (third ** fourth)) ** 0.5

if __name__ == "__main__":
    eng = Engine(ch_pressure=5, mdot=1, ch_radius=0.06)
    print(eng.results())
    for k, v in zip(eng.__dict__.keys(),eng.__dict__.values()):
        print(k, "=", v)
    eng.plot_engine()