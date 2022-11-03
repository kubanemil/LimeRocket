"""
1. Make different classes for: 
    Tanks - pressure loss in the engine working process, their thickness due to pressure, and their mass
    Engine - mdot, its mass
2. Make them the part of the class Rocket, so you can access them by:
    Rocket.Engine.mdot() or Rocket.Tanks.oxygen.mass()
"""
from engine import Engine
from tanks import Tanks
import flight
from tools import *


class Rocket:
    def __init__(self, tanks_height=1.5, radius=0.05, others_mass=5, mdot=0.5, o_f=1):
        self.gap_width = 5 * (10**(-3))
        self.engine = Engine(ch_pressure=5, mdot=mdot, ch_radius=radius-self.gap_width, o_f=o_f)
        self.tanks = Tanks(radius=radius-self.gap_width, tanks_height=tanks_height, o_f=o_f, mdot=mdot)
        self.dens = 7900  # kg/m3
        self.mdot = mdot
        self.radius = radius
        self.thick = 0.001  # meters
        self.cap_h = 0.25
        self.tanks_height = tanks_height
        self.others = others_mass  # kg
        self.engineFuel_h = 0.1
        self.fuelLox_h = 0.1
        self.gasCap_h = 0.1
        self.loxGas_h = 0.1
        self.tank_radius = radius*0.95
        self.ch_radius = radius*0.95
        self.base_h = (self.gasCap_h + self.tanks.gas_height + self.loxGas_h + self.fuelLox_h + self.tanks.total_height +
                       self.engineFuel_h + self.engine.total_length)
        self.fings_h = self.base_h/4
        self.right_h = self.fings_h/2
        self.fings_w = self.radius*3
        self.fings_num = 3
        self.max_mass = 40  # kg
        self.lox_tube = self.tanks.lox_height + ((self.loxGas_h+self.fuelLox_h)/2)
        self.fuel_tube = self.tanks.fuel_height + ((self.engineFuel_h + self.fuelLox_h) / 2)
        self.engine_tube = self.engine.total_length + (self.engineFuel_h/2)
        self.base_mass = (2*pi*self.radius * self.base_h * self.thick)*self.dens
        self.fings_mass = (self.fings_num*self.thick*self.fings_w*(self.fings_h+self.right_h)/2)*self.dens
        self.cap_mass = ((pi*self.radius * ((self.cap_h**2 + self.radius**2)**0.5))*self.thick)*self.dens
        self.check()

    def check(self):
        if self.total_mass() > self.max_mass:
            print(round(self.total_mass(), 2), "kg. ", f"Rocket heavier than {self.max_mass} kg!")
            return False
        if self.total_mass() >= self.engine.real_thrust/9.8:
            # raise ValueError("Engine can't lift the rocket!")
            print(round(self.total_mass(), 2), "kg.", round(self.engine.real_thrust/9.8, 2), "kg. ", "Engine can't lift it!")
            print(round(self.radius*100, 2), "cm.")
            return False
        return True

    def total_mass(self):
        return self.base_mass + self.cap_mass + self.fings_mass + self.engine.total_mass() + \
               self.tanks.total_mass() + self.others

    def k(self, Cd=0.75):
        air_dens = 1.2  # kg/m3
        return (Cd * air_dens * pi / 2) * (self.radius**2)

    def burn_time(self):
        return self.tanks.prop_mass()/self.mdot

    def l_r(self):
        return {
            "length": {
                "base": scaled_list([0, 0, self.base_h, self.cap_h]),
                "fings": [0, 0, self.right_h, self.fings_h],
                "engine": self.engine.reversed_lengths,
                "fuel_tank": scaled_list([self.engine.total_length+self.engineFuel_h, 0, self.tanks.fuel_height, 0]),
                "lox_tank": scaled_list([self.engine.total_length+self.engineFuel_h + self.tanks.fuel_height +
                                         self.fuelLox_h, 0, self.tanks.lox_height, 0]),
                "gas_tank": scaled_list([self.engine.total_length+self.engineFuel_h + self.tanks.fuel_height +
                                         self.fuelLox_h + self.tanks.lox_height + self.loxGas_h,
                                         0, self.tanks.gas_height, 0])
            },

            "radius": {
                "base": [0, self.radius, self.radius, 0],
                "fings": [0, self.fings_w + self.radius, self.fings_w+self.radius, self.radius],
                "engine": self.engine.reversed_radiuses,
                "tank": [0, self.tank_radius, self.tank_radius, 0],
            }
        }

    def plot_rocket(self, ix=999, save=False):
        l_r = self.l_r()
        rads = l_r['radius']
        lens = l_r['length']

        plt.rcParams["figure.figsize"] = [13.9, 6.5]
        plt.rcParams["figure.autolayout"] = True
        plt.figure(ix)
        plt.axis([0, 280, -40, 40])
        plt.axis("scaled")
        plot_part(lens['fings'], rads['fings'], fill_color="gray")
        plot_part(lens['base'], rads['base'], fill_color="black")
        plot_part(lens['engine'], rads['engine'], fill_color="red")
        plot_part(lens['fuel_tank'], rads['tank'], fill_color="darkblue")
        plot_part(lens['lox_tank'], rads['tank'], fill_color="darkred")
        plot_part(lens['gas_tank'], rads['tank'], fill_color="blue")
        plt.xlabel("Length, cm")
        plt.ylabel("Radius, cm")
        plt.figtext(0.06, 0.89, s=f"Thrust={round(self.engine.real_thrust/9.8, 2)} kg", fontsize='12')
        plt.figtext(0.06, 0.85, s=f"ISP={round(self.engine.isp, 2)} sec", fontsize='12')
        plt.figtext(0.21, 0.89, s=f"Total mass={round(self.total_mass(), 2)} kg", fontsize='12')
        plt.figtext(0.21, 0.85, s=f"Propellant mass={round(self.tanks.prop_mass(), 2)} kg", fontsize='12')
        plt.figtext(0.38, 0.89, s=f"Tanks' height={round(self.tanks_height, 2)} m", fontsize='12')
        plt.figtext(0.38, 0.85, s=f"Burn time={round(self.burn_time(), 2)} sec", fontsize='12')
        plt.figtext(0.55, 0.89, s=f"Po_Pf={round(self.tanks.Po_Pf, 2)}", fontsize='12')
        plt.figtext(0.55, 0.85, s=f"Max drag height={round(flight.Flight(self).max_height())} m", fontsize='12')
        gas_x = l_r['length']['gas_tank'][1]/2.4
        lox_x = l_r['length']['lox_tank'][1]/2.2
        fuel_x = lox_x/1.9
        print(lox_x, fuel_x)
        plt.figtext(gas_x, 0.49, s=f"Gas", fontsize='10')
        plt.figtext(lox_x, 0.49, s=f"LOX", fontsize='10')
        plt.figtext(fuel_x, 0.49, s=f"Fuel", fontsize='10')
        plt.figtext(0.07, 0.49, s=f"Engine", fontsize='9')

        if save:
            plt.savefig(f"rocket_plots\\rocket{ix}.png")
        else:
            plt.show()


if __name__ == "__main__":
    roc = Rocket(radius=0.07, tanks_height=1.2, mdot=0.5)
    for k, v in zip(roc.__dict__.keys(), roc.__dict__.values()):
        print(k, v)
    print("Base", roc.base_mass, "Cap", roc.cap_mass, "Fings", roc.fings_mass,
          "Tanks", roc.tanks.total_mass(), "Engine", roc.tanks.total_mass(), "other", roc.others)
    tank = roc.tanks
    print(tank.lox_height, tank.fuel_height, tank.tanks_mass(), tank.fuel_mass(), tank.lox_mass())
    print("TUBES:", roc.lox_tube, roc.fuel_tube, roc.engine_tube)
    print(roc)
    roc.plot_rocket()
