from math_tools import *

class Tanks:
    def __init__(self, radius=0.07):
        self.ox_height = 0.5
        self.fuel_height = 0.5
        self.pres_height = 0.5
        self.radius = radius
        self.thick = 10**(-3)  # thickness of the walls.
        self.mdot = 0.5  # kg/s

    # TODO: include pressure to mdot equation.
    # todo: returns pressure inside pressure tank.
    def lox_mass(self):
        dens = 1141     # kg/m3
        area = R_to_A(self.radius)
        h = self.ox_height
        per = 0.8  # part of tank it will take
        return per*dens*area*h

    def fuel_mass(self):
        dens = 800     # kg/m3
        area = R_to_A(self.radius)
        h = self.fuel_height
        per = 0.8  # part of tank it will take
        return per*dens*area*h

    def lox_param(self):
        return {
            "L" : 6820,  # J/mol heat of vapor for O2
            "R" : 8.31,  # J/K.mol
            "Po" : 100000,  # Pascal reference pressure
            "To" : 90,  # Kelvin reference temperature
        }

    def lox_boil_pres(self, T):
        data = self.lox_param()
        L = data["L"]; R = data["R"]; To = data["To"]; Po = data['Po']
        inx = (L / R) * ((1 / To) - (1 / T))
        return Po * (e ** inx)

    def lox_boil_temp(self, P):
        data = self.lox_param()
        L = data["L"]; R = data["R"]; To = data["To"]; Po = data['Po']
        inx = (R/L) * log(P/Po)
        return 1/ ((1/To)-inx)



tank1 = Tanks()
p = tank1.lox_boil_pres(T=120)
t = tank1.lox_boil_temp(P=10*(10**5))
print(p/(10**5))
print(t)