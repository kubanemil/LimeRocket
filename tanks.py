from tools import *

class Tanks:
    def __init__(self, radius, tanks_height, o_f, mdot, Po_Pf, ox_dens):
        self.o_f = o_f
        self.ox_dens = ox_dens
        self.total_height = tanks_height
        self.fuel_height = self.total_height/((800*self.o_f/ox_dens)+1)
        self.ox_height = (800/ox_dens)*self.o_f*self.fuel_height
        self.thick = 1 * 1e-3  # thickness of the walls.
        self.radius = radius-self.thick*3
        self.mdot = mdot  # kg/s
        self.mdot_fuel = mdot/(1+o_f)
        self.mdot_ox = mdot - self.mdot_fuel
        self.per = 0.9  # part of tank it will take

        self.pipe_radius = 4 * 1e-3  # radius of supply pipe of Fuel and OX
        self.pipe_length_fuel = 0.9   # length of supply pip of Fuel
        self.pipe_length_ox = 1.6  # length of supply pip of OX

        self.inj_tube_radius = 2 * 1e-2
        self.inj_gap = 5 * 1e-3
        self.gap_area = np.pi * self.inj_gap * (2*self.inj_tube_radius - self.inj_gap)

        self.N = 8  # number of orifices
        self.orifice_rad_fuel = 1*1e-3  # meters
        self.orifice_rad_ox = 1*1e-3  # meters
        self.inj_area_fuel = (np.pi*self.orifice_rad_fuel**2)*8  # total area of all fuel orifices
        self.inj_area_ox = (np.pi*self.orifice_rad_ox**2)*8  # total area of all oxidizer orifices
        self.inj_length_ox = 3 * 1e-3
        self.inj_length_fuel = 3 * 1e-3

        self.v_ox = self.orf_vel(self.mdot_ox, self.ox_dens, self.inj_area_ox)  # m/s velocity of OX from injector
        self.v_fuel = self.orf_vel(self.mdot_fuel, 800, self.inj_area_fuel) # m/s velocity of FUEL from injector

        self.p_diff_gap_ox = self.p_diff("L", A2R(self.gap_area), self.inj_length_ox)
        self.p_diff_pipe_fuel = self.p_diff("F", self.pipe_radius, self.pipe_length_fuel)
        self.p_diff_pipe_ox = self.p_diff("L", self.pipe_radius, self.pipe_length_ox)
        self.p_diff_inj_fuel = self.p_diff("F", A2R(self.inj_area_fuel), self.inj_length_fuel)
        self.p_diff_inj_ox = self.p_diff("L", A2R(self.inj_area_ox), self.inj_length_fuel)
        self.p_diff_fuel = self.p_diff_pipe_fuel + self.p_diff_inj_fuel
        self.p_diff_ox = self.p_diff_pipe_ox + self.p_diff_gap_ox + self.p_diff_inj_ox

        self.Po_Pf = Po_Pf
        self.gas_height = (self.total_height*(1 - self.Po_Pf*(1-self.per)))/(self.Po_Pf-1)
        self.total_mass = self.tanks_mass() + self.fuel_mass() + self.ox_mass()
        self.prop_mass = self.fuel_mass() + self.ox_mass()

    def p_diff(self, L_F_W, r, l):
        if L_F_W == "L":
            w = 0.00125  # Pa*s dynamic viscosity of hydrogen peroxide
            d = self.ox_dens  # kg/s
            mdot = self.mdot_ox
        if L_F_W == "F":
            w = 0.00164  # Pa*s dynamic viscosity of kerosene
            d = 800  # kg/s
            mdot = self.mdot_fuel
        if L_F_W == "W":
            w = 0.001  # Pa*s dynamic viscosity of water
            d = 1000  # kg/s
            mdot = self.mdot
        return (8 * w * l * mdot) / (3.14 * d * r ** 4)

    @staticmethod
    def orf_vel(mdot, dens, orf_area):
        return mdot/(dens*orf_area)

    @staticmethod
    def orifice_radius(area, N=8):
        return (area/(N*np.pi)) ** 0.5

    def ox_mass(self):
        dens = self.ox_dens     # kg/m3
        area = R2A(self.radius)
        h = self.ox_height
        return self.per*dens*area*h

    def fuel_mass(self):
        dens = 800     # kg/m3
        area = R2A(self.radius)
        h = self.fuel_height
        return self.per*dens*area*h

    def tanks_mass(self):
        perimeter = 2*pi*self.radius
        A = R2A(self.radius)
        f_v = (perimeter*self.fuel_height + 2*A) * self.thick
        ox_v = (perimeter*self.ox_height + 2*A) * self.thick
        g_v = (perimeter*self.gas_height + 2*A) * self.thick
        return 7900*(f_v+ox_v+g_v)


if __name__ == "__main__":
    pass

