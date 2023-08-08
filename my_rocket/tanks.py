from tools import *


class Tanks:
    def __init__(self, o_f, mdot, Po, tank_per=0.7, radius=0.06, tank_height=0.5):
        self.o_f = o_f
        self.Po = Po  # Initial pressure in gas and ox/fuel tanks.
        self._p_diff = 3  # bar
        self.ox_dens = 1_450_000/((38*-4.50)+1450)
        self.ox_height = tank_height
        self.fuel_height = tank_height
        self.gas_height = tank_height*2
        self.total_height = tank_height*2
        self.radius = radius
        self.mdot = mdot  # kg/s
        self.mdot_fuel = mdot/(1+o_f)
        self.mdot_ox = mdot - self.mdot_fuel
        self.ox_per = tank_per  # part of the oxidizer tank it will take
        self.fuel_per = tank_per  # part of fuel tank it will take
        self.thick = 1e-3
        self.Pf = self.Pf()
        self.viscosity_coeff = 8

        self.pipe_radius = 4 * 1e-3  # radius of supply pipe of Fuel and OX
        self.pipe_length_fuel = 0.5   # length of supply pip of Fuel
        self.pipe_length_ox = 0.5  # length of supply pip of OX

        self.N = 8  # number of orifices
        self.orifice_rad_fuel = 10*1e-4  # meters
        self.orifice_rad_ox = 10*1e-4  # meters
        self.inj_area_fuel = (np.pi*self.orifice_rad_fuel**2)*8  # total area of all fuel orifices
        self.inj_area_ox = (np.pi*self.orifice_rad_ox**2)*8  # total area of all oxidizer orifices
        self.inj_length_ox = 1*1e-3
        self.inj_length_fuel = 1*1e-3

        self.p_diff_fuel = self.p_diff("W", self.pipe_radius, self.pipe_length_fuel,
                                       A2R(self.inj_area_fuel), self.inj_length_fuel)
        self.p_diff_ox = self.p_diff("W", self.pipe_radius, self.pipe_length_ox,
                                     A2R(self.inj_area_ox), self.inj_length_fuel)

        self.mdot_fuel_total = self.p_diff_to_mdot("F", self.pipe_radius, self.pipe_length_fuel,
                                                   A2R(self.inj_area_fuel), self.inj_length_fuel,)
        self.mdot_ox_total = self.p_diff_to_mdot("L", self.pipe_radius, self.pipe_length_fuel,
                                                 A2R(self.inj_area_ox), self.inj_length_fuel)
        self.mdot_total = self.mdot_fuel_total + self.mdot_ox_total

        self.total_mass = self.tanks_mass() + self.fuel_mass() + self.ox_mass()
        self.prop_mass = self.fuel_mass() + self.ox_mass()

    def p_diff(self, L_F_W, r_pipe, l_pipe, r_inj, l_inj):  # pressure difference to create a necessary mass flow
        if L_F_W == "L":
            w = 0.00125  # Pa*s dynamic viscosity of hydrogen peroxide
            d = self.ox_dens  # kg/s
            mdot = self.mdot/2
        if L_F_W == "F":
            w = 0.00164  # Pa*s dynamic viscosity of kerosene
            d = 800  # kg/s
            mdot = self.mdot/2
        if L_F_W == "W":
            w = 0.001  # Pa*s dynamic viscosity of water
            d = 1000  # kg/s
            mdot = self.mdot
        k = self.viscosity_coeff
        l_coeff = ((l_pipe / r_pipe ** 4) + (l_inj / r_inj ** 4))
        return (k * w * mdot * l_coeff) / (pi * d)

    def p_diff_to_mdot(self, L_F_W, r_pipe, l_pipe, r_inj, l_inj):
        if L_F_W == "L":
            w = 0.00125  # Pa*s dynamic viscosity of hydrogen peroxide
            d = self.ox_dens  # kg/s
        if L_F_W == "F":
            w = 0.00164  # Pa*s dynamic viscosity of kerosene
            d = 800  # kg/s
        if L_F_W == "W":
            w = 0.001  # Pa*s dynamic viscosity of water
            d = 1000  # kg/s
        k = self.viscosity_coeff
        l_coeff = ((l_pipe / r_pipe ** 5) + (l_inj / r_inj ** 5))
        p_diff = self._p_diff * 1e5
        # print('pipe', (p_diff * pi * d)/(k * w * (l_pipe / r_pipe ** 4)))
        # print('inj', (p_diff * pi * d) / (k * w * (l_inj / r_inj ** 4)))
        return (p_diff * pi * d)/(k * w * l_coeff)

    def ox_mass(self):
        dens = self.ox_dens     # kg/m3
        area = R2A(self.radius)
        h = self.ox_height * self.ox_per  # height that filled with oxidizer
        return dens*(area*h)

    def fuel_mass(self):
        dens = 800     # kg/m3
        area = R2A(self.radius)
        h = self.fuel_height * self.fuel_per  # height that filled with fuel
        return dens*(area*h)

    def tanks_mass(self):
        perimeter = 2*pi*self.radius
        A = R2A(self.radius)
        f_v = (perimeter*self.fuel_height + (2*A)) * self.thick
        ox_v = (perimeter*self.ox_height + (2*A)) * self.thick
        g_v = (perimeter*self.gas_height + (2*A)) * self.thick
        return 7900*(f_v+ox_v+g_v)

    def Pf(self):
        A = pi * self.radius**2
        Vog = A * self.gas_height
        Vol = A * (1-self.ox_per) * self.ox_height
        Vof = A * (1-self.fuel_per) * self.fuel_height
        Vo = Vog + Vol + Vof
        Vf = A*(self.gas_height+self.ox_height+self.fuel_height)
        Pf = self.Po * (Vo/Vf)
        return Pf


if __name__ == "__main__":
    pass

