"""
The class model of rocket engine that returns CEA results, also
estimating the sizing of the engine.
"""
# TODO: Find diverging angle
# TODO: Create a method to draw a rocket scheme with valid proportions
# TODO: Estimate how fast will the rocket reach to it's critical point of heating.
# TODO: Estimate the thickness of the engine.

from CEA_Wrap import Fuel, Oxidizer, RocketProblem, DataCollector
from math_tools import *
from engine2 import Ae_At, th_area

#!!! Create code to find the area ratios!!
class Engine:
    def __init__(self, ch_pressure=5, mdot=0.5, ch_radius=0.1):
        self.ch_pressure = ch_pressure #bar
        self.mdot = mdot #kg/s
        self.o_f = 1 #oxidizer to fuel ratio
        self.ch_radius = ch_radius #m
        self.ch_area = R_to_A(self.ch_radius) #m2
        self.thrust = self.exit_vel*self.mdot

        self.create_properties()
        self.create_engine_sizes()

    def create_properties(self):
        results = self.results()
        self.exit_vel = results.son * results.mach; self.exit_son = results.son; self.exit_mach = results.mach;
        self.exit_pressure = results.p; self.exit_temp = results.t; self.exit_dens = results.rho
        self.ch_pressure = results.c_p; self.ch_temp = results.c_t; self.ch_dens = results.c_rho
        self.th_pressure = results.t_p; self.th_temp = results.t_t; self.th_dens = results.t_rho
        self.gamma = results.gamma; self.molar_mass = results.mw/1000; self.isp = results.isp

    def create_engine_sizes(self):
        self.Ae_At = Ae_At(self)
        self.th_area = th_area(self)
        self.thickness = 2*(10**(-3)) #meters
        self.exit_area = self.Ae_At*self.th_area
        self.exit_radius = A_to_R(self.exit_area)
        self.th_radius = A_to_R(self.th_area)

        self.exit_length = (self.exit_radius-self.th_radius)/(tan(15*pi/180))
        self.th_length = (self.ch_radius - self.th_radius) / (tan(30 * pi / 180))
        self.ch_length = 0.1 #meters

        self.exit_nozzle_surarea = self.nozzle_surface_area(self.exit_radius, a=15)
        self.th_nozzle_surarea = self.nozzle_surface_area(self.ch_radius, a=30)
        self.ch_surarea = 2*pi*self.ch_radius*self.ch_length

    #will return all cea results of engine, except chemical exhausts
    def results(self):
        results = RocketProblem(pressure=self.ch_pressure, materials=[Fuel('RP-1'), Oxidizer('O2(L)', temp=100)], \
                                pressure_units="bar", o_f=self.o_f, pip=self.ch_pressure).run()
        results.set_fac_ma = self.mdot / self.ch_area  # kg/s / m2

        for el in list(results.keys())[:12]:
            del results[el]

        return results


    def ch_surface_area(self):
        pass

    def nozzle_surface_area(self, R, a):
        tana = tan(a*pi/180) #degree divergence angle
        throat_area = A_to_R(self.th_area)
        return (pi/tana)*((R**2) - (throat_area**2))


    def total_mass(self):
        pass

    def __repr__(self):
        return f"Engine's Thrust: {round(self.thrust/9.8, 2)} kg.\n" + \
               f"Engine's exit velocity: {round(self.exit_vel, 2)} m/s.\n" + \
            f"Engine's chamber pressure and temperature: {self.ch_pressure} bar, {self.ch_temp} K.\n" + \
            f"Specific impulse, sea-level: {round(self.isp, 2)} seconds."

engine1 = Engine()
print(engine1.results())
print(engine1.exit_mach, engine1.gamma, engine1.exit_son)
print(engine1.Ae_At)

print(engine1.th_area ** 10000)
print(engine1.exit_nozzle_surarea)
