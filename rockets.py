from math import sqrt

import matplotlib.pyplot as plt

"""
1. Make different classes for: 
    Tanks - pressure loss in the engine working process, their thiknes due to pressure, and their mass
    Engine - mdot, its mass
2. Make them the part of the class Rocket, so you can access them by:
    Rocket.Engine.mdot() or Rocket.Tanks.oxygen.mass()
"""
from engine import Engine
from tanks import Tanks
from tools import *
from random import shuffle

class Rocket():
    def __init__(self, tanks_height=2, radius=0.07, others_mass=4, mdot=1):
        self.mdot = mdot
        self.radius = radius
        self.lthick = 0.001 #meters
        self.cap_h = 0.15;
        self.tanks_height = tanks_height
        self.others = others_mass #kg
        self.engineTank_h = 0.1; self.fuelLox_h = 0.05; self.loxCap_h = 0.05
        self.tank_radius = radius*0.95
        self.ch_radius = radius*0.95
        self.base_h = (self.loxCap_h +self.fuelLox_h + self.Tanks().total_height \
                       + self.engineTank_h + self.Engine().total_length)
        self.fing_h = self.base_h/3
        self.right_h = self.fing_h/3
        self.fing_w = self.radius*2
        self.fing_num = 4
        self.max_mass = 50 #kg
        # self.check()

    def Tanks(self):
        return Tanks(radius=self.tank_radius, tanks_height=self.tanks_height)

    def Engine(self):
        return Engine(ch_radius=self.ch_radius, mdot=self.mdot)

    def check(self):
        if self.total_mass() > self.max_mass:
            print(round(self.total_mass(), 2), "kg. ", f"Rocket heavier than {self.max_mass} kg!")
            return False
        if self.total_mass() >= self.Engine().thrust/9.8:
            # raise ValueError("Engine can't lift the rocket!")
            print(round(self.total_mass(), 2), "kg. ", round(self.Engine().thrust/9.8, 2), "kg. ", "Engine can't lift it!")
            print(round(self.radius*100, 2), "cm.")
            return False
        return True
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

    def burntime(self):
        return self.Tanks().prop_mass()/self.mdot

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

    def plot_rocket(self, ix=999, save=False):
        LRs = self.LRs()
        rads = LRs['radius']
        lens = LRs['length']
        rad_max = max([max(rads[r]) for r in rads])*100
        len_max = (self.base_h+self.cap_h)*100

        plt.rcParams["figure.figsize"] = [13.9, 6.5]
        plt.rcParams["figure.autolayout"] = True
        plt.figure(ix)
        plt.axis([0, 280, -40, 40])
        plt.axis("scaled")
        plot_part(lens['fing'], rads['fing'], fill_color="gray")
        plot_part(lens['base'], rads['base'], fill_color="black")
        plot_part(lens['engine'], rads['engine'], fill_color="red")
        plot_part(lens['fuel_tank'], rads['fuel_tank'], fill_color="darkblue")
        plot_part(lens['lox_tank'], rads['lox_tank'], fill_color="darkred")
        plt.xlabel("Length, cm")
        plt.ylabel("Radius, cm")
        plt.figtext(0.06, 0.89, s=f"Thrust={round(self.Engine().thrust/9.8, 2)} kg", fontsize='12')
        plt.figtext(0.06, 0.85, s=f"ISP={round(self.Engine().isp, 2)} sec", fontsize='12')
        plt.figtext(0.21, 0.89, s=f"Total mass={round(self.total_mass(), 2)} kg", fontsize='12')
        plt.figtext(0.21, 0.85, s=f"Propellant mass={round(self.Tanks().prop_mass(), 2)} kg", fontsize='12')
        plt.figtext(0.38, 0.89, s=f"Tanks' height={round(self.tanks_height, 2)} m", fontsize='12')
        plt.figtext(0.38, 0.85, s=f"Burntime={round(self.burntime(), 2)} sec", fontsize='12')
        lox_X = self.LRs()['length']['lox_tank'][1]/2
        fuel_X = lox_X/1.9
        print(lox_X, fuel_X)
        plt.figtext(lox_X, 0.49, s=f"LOX", fontsize='10')
        plt.figtext(fuel_X, 0.49, s=f"Fuel", fontsize='10')
        plt.figtext(0.07, 0.49, s=f"Engine", fontsize='9')

        if save:
            plt.savefig(f"rocket_plots\\rocket{ix}.png")
        else:
            plt.show()

if __name__ == "__main__":
    # hs = [h/10 for h in range(5, 19)]
    # rs = [r/100 for r in range(2, 10) ]
    # shuffle(hs)
    # rocs = []
    # for r in rs:
    #     for h in hs[:10]:
    #         rocs.append(Rocket(tanks_height=h, radius=r))
    #
    # for i in range(50):
    #     shuffle(rocs)
    #     rocs[i].plot_rocket(i, save=True)
    roc = Rocket()
    roc.plot_rocket()