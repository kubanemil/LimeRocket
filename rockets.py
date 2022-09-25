from math import sqrt
"""
1. Make different classes for: 
    Tanks - pressure loss in the engine working process, their thiknes due to pressure, and their mass
    Engine - mdot, its mass
2. Make them the part of the class Rocket, so you can access them by:
    Rocket.Engine.mdot() or Rocket.Tanks.oxygen.mass()
"""
from engine import Engine

class Rocket():
    def __init__(self, tank_h=50, ox_is_liquid=True, lox_to_fuel=2, width=14):
        self.width = width
        self.lthick = 2
        self.ox_is_liquid = ox_is_liquid
        self.cap_h = 15; self.cap_w = width
        self.tank_h = tank_h; self.tank_w = width-1
        self.engine_h = 20; self.engine_w = width-1
        self.fing_h = 35; self.fing_w = 15; self.fing_num = 3
        self.lox_to_fuel = lox_to_fuel
        self.others = 4; self.other_h = 30; self.other_w = width
        self.base_h = (self.tank_h + self.other_h + self.engine_h); self.base_w = width

    def mass(self, vol_cm, d=7.8):
        return (vol_cm / 1000) * d

    def Base(self):
        volume_cm3 = 2*3.14*(self.base_w/2) * self.base_h * (self.lthick/10)
        return self.mass(volume_cm3)

    def Cap(self):
        volume_cm3 = ((3.14*self.cap_h*self.cap_w)/2) * (self.lthick/10)
        return self.mass(volume_cm3)

    def Fings(self):
        volume_cm3 = (((self.fing_h*self.fing_w)/2)*(self.lthick/10))*self.fing_num
        return self.mass(volume_cm3)

    def Tank(self):
        volume_cm3 = 2 * 3.14 * (self.tank_w / 2) * self.tank_h * (self.lthick / 10)
        volume_cm3 += 3.14*((self.tank_w/2)**2)*(self.lthick/10)*4
        return self.mass(volume_cm3)

    def Engine(self, thick_mm=2):
        volume_cm3 = 2 * 3.14 * (self.engine_w / 2) * self.engine_h * (thick_mm / 10)
        volume_cm3 += 3.14*((self.engine_w/2)**2)*(thick_mm/10)
        return self.mass(volume_cm3)


    def lox_fuel(self):
        fuel_dens = 0.8;
        if self.ox_is_liquid:
            lox_dens = 1.14
        else:
            lox_dens = 0.0256
        inx = ((self.lox_to_fuel*fuel_dens)/lox_dens) + 1
        Area = 3.14* ((self.tank_w/2)**2)
        m_fuel = (Area * fuel_dens * self.tank_h)/inx
        m_lox = m_fuel*self.lox_to_fuel
        return [m_lox/1000, m_fuel/1000]

    def Prop(self):
        return sum(self.lox_fuel())

    def hlox_to_hfuel(self):
        fuel_dens = 0.8; lox_dens = 1.14
        return ((self.lox_to_fuel * fuel_dens) / lox_dens)

    def total(self):
        return self.Base() + self.Cap() + self.Fings() + self.Engine()\
                +self.Tank() + self.Prop() + self.others

    def K(self, drag_const=0.75):
        return 0.5 * 1.2 * 3.14 * ((self.width/2) ** 2) * drag_const * (10**(-4))



if __name__ == "__main__":
    rocket = Rocket(tank_h=100, ox_is_liquid=False)
    print(rocket.K())