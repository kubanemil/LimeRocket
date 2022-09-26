"""
The class model of rocket engine that returns CEA results, also
estimating the sizing of the engine.
"""
# TODO: Create a method to draw a rocket scheme with valid proportions
# TODO: Estimate how fast will the rocket reach to it's critical point of heating.

from math_tools import *
from CEA_Wrap import Fuel, Oxidizer, RocketProblem
from engine2 import create_properties, create_engine_sizes

class Engine:
    def __init__(self, ch_pressure=5, mdot=0.5, ch_radius=0.1):
        self.ch_pressure = ch_pressure  # bar
        self.mdot = mdot    # kg/s
        self.o_f = 1    # oxidizer to fuel ratio
        self.ch_radius = ch_radius  # m
        self.ch_area = R_to_A(self.ch_radius)   # m2
        self.ch_length = 0.1
        create_properties(self)
        create_engine_sizes(self)
        self.thrust = self.exit_vel * self.mdot

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

    def __repr__(self):
        return f"Thrust: {round(self.thrust/9.8, 2)} kg.\n" + \
               f"Exit velocity: {round(self.exit_vel, 2)} m/s.\n" + \
            f"Chamber pressure and temperature: {self.ch_pressure} bar, {self.ch_temp} K.\n" + \
            f"Specific impulse, sea-level: {round(self.isp, 2)} seconds."


