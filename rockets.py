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
    def __init__(self, tanks_height=1, radius=0.07, others=4):
        self.radius = radius
        self.lthick = 0.001 #meters
        self.cap_h = 0.15;
        self.fing_h = 0.45; self.right_h=0.25; self.fing_w = 0.15; self.fing_num = 3
        self.tanks_height = tanks_height
        self.others = others #kg
        self.engineTank_h = 0.1; self.fuelLox_h = 0.05; self.loxCap_h = 0.05
        self.tank_radius = 0.065
        self.ch_radius = 0.06
        self.base_h = (self.loxCap_h +self.fuelLox_h + self.Tanks().total_height \
                       + self.engineTank_h + self.Engine().total_length)
        self.check_thrust()

    def Tanks(self):
        return Tanks(radius=self.tank_radius, tanks_height=self.tanks_height)

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
        volume = ((((self.right_h*self.fing_w))+((self.fing_h-self.right_h)*self.fing_w/2))*(self.lthick))*self.fing_num
        return self.mass(volume)

    def total_mass(self):
        return self.Base() + self.Cap() + self.Fings() + self.Engine().total_mass()\
                +self.Tanks().total_mass() + self.others

    def K(self, drag_const=0.8):
        return 0.5 * 1.2 * 3.14 * (self.radius ** 2) * drag_const * (10**(-4))

    def all_parts(self):
        return [self.Engine().total_mass(), self.Tanks().total_mass(), self.Base(), self.Cap(), self.Fings(), self.others]

    def LRs(self):
        return {
            "length":{
                "base" : scaled_list([0, 0, self.base_h, self.cap_h]),
                "fing": [0, 0, self.right_h, self.fing_h],
                "engine": self.Engine().reversed_lengths,
                "fuel_tank" : scaled_list([self.Engine().total_length+self.engineTank_h, 0, self.Tanks().fuel_height, 0]),
                "lox_tank" : scaled_list([self.Engine().total_length+self.engineTank_h+
                                           self.Tanks().fuel_height+self.fuelLox_h,0, self.Tanks().lox_height, 0]),
            },

            "radius" : {
                "base" : [0, self.radius, self.radius, 0],
                "fing": [0, self.fing_w + self.radius, self.fing_w+self.radius, self.radius],
                "engine": self.Engine().reversed_radiuses,
                "fuel_tank": [0, self.tank_radius, self.tank_radius, 0],
                "lox_tank": [0, self.tank_radius, self.tank_radius, 0],
            }
        }


    def plot_rocket(self):
        LRs = self.LRs()
        rads = LRs['radius']
        lens = LRs['length']
        rad_max = max([max(rads[r]) for r in rads])*100
        len_max = (self.base_h+self.cap_h)*100

        plt.axis([0, len_max+10, -rad_max-10, rad_max+10])
        plt.axis("scaled")
        plot_part(lens['fing'], rads['fing'], fill_color="white")
        plot_part(lens['base'], rads['base'], fill_color="white")
        plot_part(lens['engine'], rads['engine'], fill_color="white")
        plot_part(lens['fuel_tank'], rads['fuel_tank'], fill_color="white")
        plot_part(lens['lox_tank'], rads['lox_tank'], fill_color="white")
        plt.xlabel("Length, meters")
        plt.ylabel("Radius, meters")
        plt.show()

if __name__ == "__main__":
    rocket = Rocket(tanks_height=1)
    print(rocket.K())
    print(rocket.total_mass())
    print(rocket.all_parts())
    print(rocket.plot_rocket())