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

        self.pipe_radius = 3 * (10**(-3))  # radius of supply pipe of Fuel and LOX
        self.pipe_length_fuel = 0.9   # length of supply pip of Fuel
        self.pipe_length_lox = 1.6  # length of supply pip of LOX

        self.injector_a = 50  # degree angle of injection
        self.v_lox = 20   # m/s velocity of LOX from injector
        self.v_fuel = self.v_lox * np.tan(self.injector_a * np.pi/180) # m/s velocity of FUEL from injector

        self.inj_radius = 2 * 10**(-2)
        self.inj_gap = 1 * 10**(-3)
        self.gap_area = np.pi * self.inj_gap * (2*self.inj_radius - self.inj_gap)

        self.inj_area_fuel = self.area_inj(vel=self.v_fuel, dens=800)
        self.inj_area_lox = self.area_inj(vel=self.v_lox, dens=1141)
        self.inj_length_lox = 3 * 10**(-2)
        self.inj_length_fuel = 3 * 10**(-3)
        self.orifice_rad_fuel = self.orifice_radius(self.inj_area_fuel)
        self.orifice_rad_lox = self.orifice_radius(self.inj_area_lox)

        self.p_diff_pipe_fuel = self.p_diff("F", self.pipe_radius, self.pipe_length_fuel)
        self.p_diff_inj_fuel = self.p_diff("F", A_to_R(self.inj_area_fuel), self.inj_length_fuel)
        self.p_diff_pipe_lox = self.p_diff("L", self.pipe_radius, self.pipe_length_lox)
        self.p_diff_gap_lox = self.p_diff("L", A_to_R(self.gap_area), self.inj_length_lox)
        self.p_diff_inj_lox = self.p_diff("L", A_to_R(self.inj_area_lox), self.inj_length_fuel)
        self.p_diff_fuel = self.p_diff_pipe_fuel + self.p_diff_inj_fuel
        self.p_diff_lox = self.p_diff_pipe_lox + self.p_diff_gap_lox + self.p_diff_inj_lox

        self.Po_Pf = Po_Pf
        self.gas_height = (self.total_height*(1 - self.Po_Pf*(1-self.per)))/(self.Po_Pf-1)
        self.total_mass = self.tanks_mass() + self.fuel_mass() + self.lox_mass()
        self.prop_mass = self.fuel_mass() + self.lox_mass()


    def p_diff(self, L_F_W, r, l, mdot=0.5):
        if L_F_W == "L":
            w = 7.71 * 10 ** (-6)  # Pa*s dynamic viscosity of lox
            d = 1141  # kg/s
        if L_F_W == "F":
            w = 0.00164  # Pa*s dynamic viscosity of kerosene
            d = 800  # kg/s
        if L_F_W == "W":
            w = 0.001  # Pa*s dynamic viscosity of water
            d = 1000  # kg/s
        return (8 * w * l * mdot) / (3.14 * d * r ** 4)

    def area_from_p_diff(self, L_F_W, p_diff):
        if L_F_W == "L":
            w = 7.71 * 10 ** (-6)  # Pa*s dynamic viscosity of lox
            d = 1141  # kg/s
        if L_F_W == "F":
            w = 0.00164  # Pa*s dynamic viscosity of kerosene
            d = 800  # kg/s
        if L_F_W == "W":
            w = 0.001  # Pa*s dynamic viscosity of water
            d = 1000  # kg/s

    @staticmethod
    def area_inj(vel, dens, mdot=0.5):
        return mdot/(dens*vel)

    @staticmethod
    def orifice_radius(area, N=16):
        return (area/(N*np.pi)) ** 0.5

    def inj_gap(self):
        pass

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

