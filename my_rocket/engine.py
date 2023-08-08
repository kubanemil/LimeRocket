"""
The class model of rocket engine that returns CEA results, also
estimating the sizing of the engine.
"""

from tools import *
from CEA_Wrap import Fuel, Oxidizer, RocketProblem
import matplotlib.pyplot as plt
import numpy as np

# CEA verified. #Sizing verified.


class Engine:
    def __init__(self, ch_pressure, mdot, o_f, KPD, ho2_per):
        self.ch_pressure = ch_pressure  # bar pressure in the engine's chamber
        self.mdot = mdot    # kg/s  propellant's mass flow
        self.o_f = o_f    # oxidizer to fuel ratio
        self.ho2_per = ho2_per
        self.thickness = 0.001
        self.ch_radius = 63 * 1e-3
        self.th_radius = 90/pi * 1e-3
        self.ex_radius = 100/pi * 1e-3
        self.ch_area = R2A(self.ch_radius)   # m2
        self.th_area = R2A(self.th_radius)  # m2
        self.ex_area = R2A(self.ex_radius)  # m2
        self.Ae_At = self.ex_area / self.th_area
        self.ch_length = 0.2  # length of the chamber
        self.angle_exit = 5   # degree
        self.angle_throat = 30    # degree
        self.angle_injector = 50
        self.insulation_thick = 0.2  # cm
        self.KPD = KPD
        self.create_properties()
        self.create_engine_sizes()

    def results(self):
        results = RocketProblem(pressure=self.ch_pressure, materials=[Fuel('RP-1', wt_percent=90),
                                                                      Oxidizer('H2O2(L)', wt_percent=self.ho2_per),
                                                                      Oxidizer('H2O', wt_percent=(100-self.ho2_per), temp=300)],
                                pressure_units="bar", o_f=self.o_f, ae_at=self.Ae_At).run()
        results.set_fac_ma = self.mdot / self.ch_area  # kg/s / m2
        for el in list(results.keys())[:12]:
            del results[el]
        return results

    def total_mass(self):
        total_surarea = self.ch_surarea + self.th_nozzle_surarea + self.exit_nozzle_surarea
        dens = 7900     # kg/m3
        t = self.thickness  # meters
        return round(total_surarea*t*dens, 3)

    def injector_line(self):
        x = (self.ch_radius * 100 - self.insulation_thick) / np.tan(self.angle_injector * np.pi / 180)
        return np.array([0, x])

    def create_properties(self):
        results = self.results()
        self.exit_vel = results.son * results.mach
        self.exit_son = results.son
        self.exit_mach = results.mach
        self.exit_pressure = results.p
        self.exit_temp = results.t
        self.exit_dens = results.rho
        self.ch_pressure = results.c_p
        self.ch_temp = results.c_t
        self.ch_dens = results.c_rho
        self.th_pressure = results.t_p
        self.th_temp = results.t_t
        self.th_dens = results.t_rho
        self.gamma = results.gamma
        self.molar_mass = results.mw / 1000
        self.isp = results.isp
        self.ideal_thrust = self.exit_vel * self.mdot
        self.thrust = self.KPD*self.ideal_thrust

    def nozzle_surface_area(self, th_radius, R, angle):
        tana = tan(angle * pi / 180)  # degree divergence angle
        return (pi / tana) * ((R ** 2) - (th_radius ** 2))

    def create_engine_sizes(self):
        self.exit_length = (self.ex_radius - self.th_radius) / (tan(self.angle_exit * pi / 180))
        self.th_length = (self.ch_radius - self.th_radius) / (tan(self.angle_throat * pi / 180))

        self.exit_nozzle_surarea = self.nozzle_surface_area(self.th_radius, self.ex_radius, angle=self.angle_exit)
        self.th_nozzle_surarea = self.nozzle_surface_area(self.th_radius, self.ch_radius, angle=self.angle_throat)
        self.ch_surarea = (2 * pi * self.ch_radius * self.ch_length) + (pi * (self.ch_radius ** 2))

        self.radiuses = np.array([self.ch_radius, self.ch_radius, self.th_radius, self.ex_radius, 0])
        self.lengths = np.array(scaled_list([0, self.ch_length, self.th_length, self.exit_length, 0]))
        self.reversed_radiuses = np.array([self.ex_radius, self.th_radius, self.ch_radius, self.ch_radius, 0])
        self.reversed_lengths = np.array(scaled_list([0, self.exit_length, self.th_length, self.ch_length, 0]))
        self.total_length = self.ch_length + self.th_length + self.exit_length

    def plot_engine(self, ix=999, save=False):
        lengths_cm = self.lengths*100
        radiuses_cm = self.radiuses*100
        plt.rcParams["figure.figsize"] = [13.9, 6.5]
        plt.rcParams["figure.autolayout"] = True
        plt.figure(ix)
        plt.axis([0, 30, -8, 8])
        plt.axis('scaled')
        plt.plot(lengths_cm, radiuses_cm + self.insulation_thick, "black")
        plt.plot(lengths_cm, -radiuses_cm - self.insulation_thick, "black")
        plt.plot(lengths_cm, radiuses_cm, "red")
        plt.plot(lengths_cm, -radiuses_cm, "red")
        plt.grid(True)
        plt.figtext(0.13, 0.93, s=f"Thrust={round(self.thrust/9.8, 2)} kg", fontsize='12')
        plt.figtext(0.13, 0.9, s=f"Engine's mass={round(self.total_mass(), 2)} kg", fontsize='12')
        plt.figtext(0.30, 0.93, s=f"Mass flow={round(self.mdot, 2)} kg/s", fontsize='12')
        plt.figtext(0.30, 0.9, s=f"O/F ratio={round(self.o_f, 2)}", fontsize='12')
        plt.figtext(0.45, 0.93, s=f"ISP={round(self.isp, 2)} sec",  fontsize='12')
        plt.figtext(0.45, 0.9, s=f"Exit velocity={round(self.exit_vel, 2)} m/s", fontsize='12')
        plt.figtext(0.7, 0.7, s=self.__repr__(), fontsize='10')

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
        return f"Thrust: {round(self.thrust/9.8, 2)} kg." "\n" \
               f"Exit velocity: {round(self.exit_vel, 2)} m/s." "\n" \
            f"Chamber pressure and temperature: {self.ch_pressure} bar, {self.ch_temp} K." "\n" \
            f"Specific impulse, sea-level: {round(self.isp, 2)} seconds." "\n" \
            f"Mass flow rate: {rnd(self.mdot)} kg/s" "\n" \
            f"Chamber area: {rnd(self.ch_area)} m2" "\n" \
            f"Throat area: {rnd(self.th_area)} m2" "\n" \
            f"Exit area: {rnd(self.ex_area)} m2" "\n" \
