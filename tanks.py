from tools import *

class Tanks:
    def __init__(self, radius=0.07, tanks_height=1, o_f=1, mdot=0.5):
        self.o_f = o_f
        self.total_height = tanks_height
        self.pres_height = 0.5
        self.radius = radius
        self.thick = 10**(-3)  # thickness of the walls.
        self.mdot = mdot  # kg/s
        self.per = 0.9  # part of tank it will take
        self.pipe_radius = 2 * (10**(-3))
        self.pipe_length = 1
        self.gas_height = 0.6
        self.Po_Pf = (self.total_height + self.gas_height)/(self.gas_height + ((1-self.per)*self.total_height))

    # TODO: include pressure to mdot equation.
    # todo: returns pressure inside pressure tank.

    def p_diff(self):
        # water: w = 0.0010518 Pas*s
        w = 0.00164     # Pa*s
        d = 800  # kg/s
        L = self.pipe_length
        R = self.pipe_radius
        mdot = self.mdot
        return (8 * w * L * mdot) / (3.14 * d * (R ** 4))

    def get_mdot(self, p_diff=5.0):  # mass rate flow
        w = 0.00164     # Pa*s
        d = 800  # kg/s
        L = self.pipe_length
        R = self.pipe_radius
        p_diff = p_diff * (10**5)
        return (3.14 * d * (R ** 4) * (p_diff)) / (8 * w * L)

    def P(self, h_diff, Po_gas, dens, ho, h_gas):
        if h_diff > ho:
            print("h_diff can't be greater than ho!")
        else:
            return Po_gas*(h_gas/(h_gas+h_diff)) + dens*10*(ho-h_diff)

    @property
    def fuel_height(self):
        return self.total_height/((800*self.o_f/1141)+1)

    @property
    def lox_height(self):
        return (800/1141)*self.o_f*self.fuel_height

    def lox_mass(self):
        dens = 1141     # kg/m3
        area = R_to_A(self.radius)
        h = self.lox_height
        return self.per*dens*area*h

    def fuel_mass(self):
        dens = 800     # kg/m3
        area = R_to_A(self.radius)
        h = self.fuel_height
        return self.per*dens*area*h

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

    def tanks_mass(self):
        perimeter = 2*pi*self.radius
        A = R_to_A(self.radius)
        f_v = (perimeter*self.fuel_height + 2*A) * self.thick
        ox_v = (perimeter*self.lox_height + 2*A) * self.thick
        g_v = (perimeter*self.gas_height + 2*A) * self.thick
        return 7900*(f_v+ox_v+g_v)

    def prop_mass(self):
        return self.fuel_mass() + self.lox_mass()

    def total_mass(self):
        return self.tanks_mass() + self.fuel_mass() + self.lox_mass()


if __name__ == "__main__":
    tank1 = Tanks()
    p = tank1.lox_boil_pres(T=120)
    t = tank1.lox_boil_temp(P=10*(10**5))

    print(tank1.P(h_diff=0, Po_gas=500000, h_gas=0.2, ho=0.5, dens=1000))
    print(tank1.P(h_diff=0.5, Po_gas=500000, h_gas=0.2, ho=0.5, dens=1000))
    print(tank1.get_mdot(p_diff=5))
    print(tank1.get_mdot(p_diff=5.01), tank1.get_mdot(p_diff=5.05))
    print(tank1.Po_Pf)
    print(tank1.total_mass())