from math import sqrt
"""
1. Make different classes for: 
    Tanks - pressure loss in the engine working process, their thiknes due to pressure, and their mass
    Engine - mdot, its mass
2. Make them the part of the class Rocket, so you can access them by:
    Rocket.Engine.mdot() or Rocket.Tanks.oxygen.mass()
"""
from engine import Engine
from tanks import Tanks
import matplotlib.pyplot as plt
from tools import *
import numpy as np

class Rocket():
    def __init__(self, fuel_height=0.5, radius=0.07, others=4):
        self.radius = radius
        self.lthick = 0.001 #meters
        self.cap_h = 0.15;
        self.fing_h = 0.45; self.fing_w = 0.25; self.fing_num = 3
        self.fuel_height = fuel_height
        self.others = others #kg
        self.engineTank_h = 0.1; self.fuelLox_h = 0.05; self.loxCap_h = 0.05;
        self.tank_radius = 0.065
        self.ch_radius = 0.06
        self.base_h = (self.loxCap_h + self.Tanks().lox_height + self.fuelLox_h + self.Tanks().fuel_height + \
                       self.engineTank_h + self.Engine().total_length)
        self.check_thrust()

    def Tanks(self):
        return Tanks(radius=self.tank_radius, fuel_height=self.fuel_height)

    def Engine(self):
        return Engine(ch_radius=self.ch_radius)

    def check_thrust(self):
        if self.total_mass() >= self.Engine().thrust/9.8:
            raise ValueError("Engine can't lift the rocket!")

    def mass(self, V, d=7900):
        return V * d

    def Base(self):
        volume = 2*3.14*self.radius * self.base_h * (self.lthick)
        return self.mass(volume)

    def Cap(self):
        volume = (3.14*self.cap_h*self.radius) * (self.lthick)
        return self.mass(volume)

    def Fings(self):
        volume = (((self.fing_h*self.fing_w)/2)*(self.lthick))*self.fing_num
        return self.mass(volume)

    def total_mass(self):
        return self.Base() + self.Cap() + self.Fings() + self.Engine().total_mass()\
                +self.Tanks().total_mass() + self.others

    def K(self, drag_const=0.8):
        return 0.5 * 1.2 * 3.14 * (self.radius ** 2) * drag_const * (10**(-4))

    def all_parts(self):
        return [self.Engine().total_mass(), self.Tanks().total_mass(), self.Base(), self.Cap(), self.Fings(), self.others]

    def lengthRadiuses(self):
        return {
            "lengths":{
                "fing": [0, 0, self.fing_h],
                "engine": self.Engine().reversed_lengths,
                "fuel_tank" : scaled_list([self.Engine().total_length+self.engineTank_h, 0, self.Tanks().fuel_height, 0]),
                "lox_tank" : scaled_list([self.Engine().total_length+self.engineTank_h+
                                           self.Tanks().fuel_height+self.fuelLox_h,0, self.Tanks().lox_height, 0]),
            },

            "radiuses" : {
                "fing": [0, self.fing_w + self.radius, self.radius],
                "engine": self.Engine().reversed_radiuses,
                "fuel_tank": [0, self.tank_radius, self.tank_radius, 0],
                "lox_tank": [0, self.tank_radius, self.tank_radius, 0],
            }

        }
    def plot_base(self):
        lengths_cm = np.array(scaled_list([0, 0, self.base_h, self.cap_h])) * 100
        radiuses_cm = np.array([0, self.radius, self.radius, 0]) * 100
        plt.plot(lengths_cm, radiuses_cm, 'black')
        plt.plot(lengths_cm, -radiuses_cm, 'black')
        plt.fill_between(lengths_cm, radiuses_cm, color="black")
        plt.fill_between(lengths_cm, -radiuses_cm, color="black")

    def plot_fings(self):
        lengths_cm = np.array([0, 0, self.fing_h])*100
        radiuses_cm = np.array([0, self.fing_w+self.radius, self.radius ])*100
        plt.plot(lengths_cm, radiuses_cm, 'black')
        plt.plot(lengths_cm, -radiuses_cm, 'black')
        plt.fill_between(lengths_cm, radiuses_cm, color="gray")
        plt.fill_between(lengths_cm, -radiuses_cm, color="gray")

    def plot_engine(self):
        lengths_cm = self.Engine().reversed_lengths * 100
        radiuses_cm = self.Engine().reversed_radiuses * 100
        plt.plot(lengths_cm, radiuses_cm, "red")
        plt.plot(lengths_cm, -radiuses_cm, "red")
        plt.plot(lengths_cm, radiuses_cm + self.Engine().insul_thick, "black")
        plt.plot(lengths_cm, -(radiuses_cm + self.Engine().insul_thick), "black")
        plt.fill_between(lengths_cm, radiuses_cm, color="red")
        plt.fill_between(lengths_cm, -radiuses_cm, color="red")

    def plot_fuel_tank(self):
        lengths_cm = np.array(scaled_list([self.Engine().total_length+self.engineTank_h, 0, self.Tanks().fuel_height, 0]))*100
        radiuses_cm = np.array([0, self.tank_radius, self.tank_radius, 0])*100
        plt.plot(lengths_cm, radiuses_cm, "darkblue")
        plt.plot(lengths_cm, -radiuses_cm, "darkblue")
        plt.fill_between(lengths_cm, radiuses_cm, color="blue")
        plt.fill_between(lengths_cm, -radiuses_cm, color="blue")

    def plot_lox_tank(self):
        tank = self.Tanks()
        lengths_cm = np.array(scaled_list([self.Engine().total_length+self.engineTank_h+
                                           self.Tanks().fuel_height+self.fuelLox_h,0, self.Tanks().lox_height, 0]))*100
        radiuses_cm = np.array([0, self.tank_radius, self.tank_radius, 0]) * 100
        plt.plot(lengths_cm, radiuses_cm, "darkred")
        plt.plot(lengths_cm, -radiuses_cm, "darkred")
        plt.fill_between(lengths_cm, radiuses_cm, color="darkred")
        plt.fill_between(lengths_cm, -radiuses_cm, color="darkred")

    def plot_rocket(self):
        plt.axis([0, (self.base_h*100+self.cap_h*100+20), -(self.radius*100+self.fing_w*100 + 10), \
                  (self.radius*100+self.fing_w*100 + 10)])
        plt.axis("scaled")
        plot_part()
        self.plot_base()
        self.plot_engine()
        self.plot_fuel_tank()
        self.plot_lox_tank()
        plt.xlabel("Length, meters")
        plt.ylabel("Radius, meters")
        plt.show()

if __name__ == "__main__":
    rocket = Rocket(fuel_height=0.6)
    print(rocket.K())
    print(rocket.total_mass())
    print(rocket.all_parts())
    print(rocket.plot_rocket())