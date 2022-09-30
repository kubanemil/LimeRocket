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

class Rocket():
    def __init__(self, fuel_height=0.5, width=0.14, others=4):
        self.width = width
        self.lthick = 0.001 #meters
        self.cap_h = 0.15;
        self.fing_h = 0.35; self.fing_w = 0.15; self.fing_num = 3
        self.fuel_height = fuel_height
        self.others = others #kg
        self.other_h = 0.30
        self.base_h = (self.Tanks().lox_height + self.Tanks().fuel_height + \
                       self.other_h + self.Engine().total_length)
        self.check_thrust()

    def check_thrust(self):
        if self.total_mass() >= self.Engine().thrust/9.8:
            raise ValueError("Engine can't lift the rocket!")

    def mass(self, V, d=7900):
        return V * d

    def Base(self):
        volume = 2*3.14*(self.width/2) * self.base_h * (self.lthick)
        return self.mass(volume)

    def Cap(self):
        volume = ((3.14*self.cap_h*self.width)/2) * (self.lthick)
        return self.mass(volume)

    def Fings(self):
        volume = (((self.fing_h*self.fing_w)/2)*(self.lthick))*self.fing_num
        return self.mass(volume)

    def Tanks(self):
        return Tanks(radius=(self.width-0.01)/2, fuel_height=self.fuel_height)

    def Engine(self):
        return Engine(ch_radius=(self.width-0.02)/2)

    def total_mass(self):
        return self.Base() + self.Cap() + self.Fings() + self.Engine().total_mass()\
                +self.Tanks().total_mass() + self.others

    def K(self, drag_const=0.8):
        return 0.5 * 1.2 * 3.14 * ((self.width/2) ** 2) * drag_const * (10**(-4))

    def all_parts(self):
        return [self.Engine().total_mass(), self.Tanks().total_mass(), self.Base(), self.Cap(), self.Fings(), self.others]

if __name__ == "__main__":
    rocket = Rocket()
    print(rocket.K())
    print(rocket.total_mass())
    print(rocket.all_parts())