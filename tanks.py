from tools import *

class Tanks:
    def __init__(self, radius, tanks_height, o_f, mdot, Po_Pf):
        self.o_f = o_f
        self.total_height = tanks_height
        self.fuel_height = self.total_height/((800*self.o_f/1141)+1)
        self.lox_height = (800/1141)*self.o_f*self.fuel_height
        self.pres_height = 0.5
        self.radius = radius
        self.thick = 10**(-3)  # thickness of the walls.
        self.mdot = mdot  # kg/s
        self.per = 0.9  # part of tank it will take
        self.pipe_radius = 3 * (10**(-3))
        self.pipe_length_fuel = 0.9
        self.pipe_length_lox = 1.6
        self.inj_radius_fuel = np.pi*0.001
        self.inj_radius_lox = np.pi*0.001
        self.inj_area_fuel = np.pi * self.inj_radius_fuel**2
        self.inj_area_lox = np.pi * self.inj_radius_lox**2
        self.inj_length = 2 * 10**(-2)
        self.Po_Pf = Po_Pf
        self.gas_height = (self.total_height*(1 - self.Po_Pf*(1-self.per)))/(self.Po_Pf-1)
        self.total_mass = self.tanks_mass() + self.fuel_mass() + self.lox_mass()
        self.prop_mass = self.fuel_mass() + self.lox_mass()

    def fuel_p_diff(self):
        w = 0.00164  # Pa*s dynamic viscosity of kerosene
        d = 800  # kg/s
        r = self.pipe_radius
        R = self.inj_radius_fuel
        l = self.pipe_length_fuel
        L = self.inj_length
        mdot = self.mdot/2
        pipe_p_diff = (8 * w * l * mdot) / (3.14 * d * r ** 4)
        injector_p_diff = (8 * w * L * mdot) / (3.14 * d * R ** 4)
        return pipe_p_diff,  injector_p_diff

    def lox_p_diff(self):
        w = 7.71 * 10 ** (-6)  # Pa*s dynamic viscosity of lox
        d = 1141  # kg/s
        r = self.pipe_radius
        R = self.inj_radius_lox
        l = self.pipe_length_lox
        L = self.inj_length
        mdot = self.mdot/2
        pipe_p_diff = (8 * w * l * mdot) / (3.14 * d * r ** 4)
        injector_p_diff = (8 * w * L * mdot) / (3.14 * d * R ** 4)
        return pipe_p_diff, injector_p_diff

    def fuel_inj_v(self):
        mdot = self.mdot/2
        return mdot/(800*self.inj_area)

    def lox_inj_v(self):
        mdot = self.mdot/2
        return mdot/(1141*self.inj_area)

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


if __name__ == "__main__":
    pass

