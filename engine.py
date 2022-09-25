"""
The class model of rocket engine that returns CEA results, also
estimating the sizing of the engine.
"""
"""
TODO:
1.Find diverging angle
2.Create a method to draw a rocket scheme with valid proportions
3.Estimate how fast will the rocket reach to it's critical point of heating.
4.Estimate the thikness of the rocket.
"""
from CEA_Wrap import Fuel, Oxidizer, RocketProblem, DataCollector
import math_tools


#!!! Create code to find the area ratios!!
class Engine:
    def __init__(self):
        self.ch_pressure = 5 #bar
        self.mdot = 0.5 #kg/s
        self.o_f = 1 #oxidizer to fuel ratio
        self.ch_radius = 0.1 #m
        self.ch_area = math_tools.R_to_A(self.ch_radius) #m2
        self.fuel = Fuel('RP-1')
        self.ox = Oxidizer('O2(L)', temp=100)
        self.create_properties()

    def create_properties(self):
        results = self.results()
        self.exit_vel = results.son * results.mach; self.exit_son = results.son; self.exit_mach = results.mach;
        self.exit_pressure = results.p; self.exit_temp = results.t; self.exit_dens = results.rho

        self.ch_pressure = results.c_p; self.ch_temp = results.c_t; self.ch_dens = results.c_rho

        self.th_pressure = results.t_p; self.th_temp = results.t_t; self.th_dens = results.t_rho

        self.gamma = results.gamma; self.molar_mass = results.mw/1000; self.isp = results.isp


    #will return all cea results of engine, except chemical exhausts
    def results(self):
        results = RocketProblem(pressure=self.ch_pressure, materials=[self.fuel, self.ox], \
                                pressure_units="bar", o_f=self.o_f, pip=self.ch_pressure).run()
        results.set_fac_ma = self.mdot / self.ch_area  # kg/s / m2

        for el in list(results.keys())[:12]:
            del results[el]

        return results



    # #for ideal gas:
    # def At_ideal(self):
    #     R = 8.314; g = 9.8
    #     Pt = self.throat_pressure*(10**5)
    #     sqrt_ix = (R*self.throat_temp)/(self.gamma*g)
    #     return (self.mdot/Pt) * (sqrt_ix**(1/2))
    #
    # def At2(self):
    #     Pt = self.throat_pressure
    #     mdot = self.mdot
    #     q = self.throat_dens
    #     y = self.gamma
    #     Cd = 0.8
    #     y_ix = (2/(y+1))**((y+1)/(y-1))
    #     sqrt_ix = (2*q*Pt*y_ix)**0.5
    #     return mdot/(Cd*sqrt_ix)

    @property
    def throat_area(self):
        Pc = self.ch_pressure * 100000; y = self.gamma
        Tc = self.ch_temp; mdot = self.mdot;
        R = 8.314/self.molar_mass

        first = mdot / Pc
        second = (R * Tc) / y
        third = (2 / (y + 1))
        fourth = (y + 1) / (y - 1)
        return first * (second / (third ** fourth)) ** 0.5


    #ratio of exit area (Ae) to throat area(At) = Ae/At
    # def ideal_Ae_At(self):
    #     Ma = self.exit_mach
    #     y = self.gamma
    #     first = 1 + (((y-1)/2)*(Ma**2))
    #     y_ix = (y+1)/2
    #     return (1/Ma)*((first/y_ix)**(y_ix/(y-1)))

    def Ae_At(self):
        k = self.gamma
        pip = 1/self.ch_pressure
        first_ix = ((k + 1) / 2) ** (1 / (k - 1)) * pip ** (1 / k)
        sqrt_ix = (k + 1) / (k - 1) * (1 - pip ** ((k - 1) / k))
        At_Ae = first_ix * (sqrt_ix) ** 0.5
        return 1/At_Ae

    @property
    def exit_area(self):
        return self.Ae_At()*self.At()

    @property
    def thrust(self):
        return self.exit_vel*self.mdot
    def __repr__(self):
        return f"Engine's Thrust: {round(self.thrust/9.8, 2)} kg.\n" + \
               f"Engine's exit velocity: {round(self.exit_vel, 2)} m/s.\n" + \
            f"Engine's chamber pressure and temperature: {self.ch_pressure} bar, {self.ch_temp} K.\n" + \
            f"Specific impulse, sea-level: {round(self.isp, 2)} seconds."

engine1 = Engine()
print(engine1.results())
print(engine1.exit_mach, engine1.gamma, engine1.exit_son)
print(engine1.Ae_At())

print(engine1.throat_area ** 10000)
print(engine1)
