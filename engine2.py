from math_tools import *


def create_properties(self):
    results = self.results()
    self.exit_vel = results.son * results.mach;
    self.exit_son = results.son;
    self.exit_mach = results.mach;
    self.exit_pressure = results.p;
    self.exit_temp = results.t;
    self.exit_dens = results.rho
    self.ch_pressure = results.c_p;
    self.ch_temp = results.c_t;
    self.ch_dens = results.c_rho
    self.th_pressure = results.t_p;
    self.th_temp = results.t_t;
    self.th_dens = results.t_rho
    self.gamma = results.gamma;
    self.molar_mass = results.mw / 1000;
    self.isp = results.isp

def nozzle_surface_area(th_radius, R, a):
    tana = tan(a*pi/180) #degree divergence angle
    return (pi/tana)*((R**2) - (th_radius**2))

def create_engine_sizes(self):
    self.Ae_At = Ae_At(self)
    self.th_area = th_area(self)
    self.thickness = 2*(10**(-3)) #meters
    self.exit_area = self.Ae_At*self.th_area
    self.exit_radius = A_to_R(self.exit_area)
    self.th_radius = A_to_R(self.th_area)

    self.exit_length = (self.exit_radius-self.th_radius)/(tan(15*pi/180))
    self.th_length = (self.ch_radius - self.th_radius) / (tan(30 * pi / 180))

    self.exit_nozzle_surarea = nozzle_surface_area(self.th_radius, self.exit_radius, a=15)
    self.th_nozzle_surarea = nozzle_surface_area(self.th_radius, self.ch_radius, a=30)
    self.ch_surarea = (2*pi*self.ch_radius*self.ch_length) + (pi*(self.ch_radius**2))


def Ae_At(self):
    k = self.gamma
    pip = 1 / self.ch_pressure
    first_ix = ((k + 1) / 2) ** (1 / (k - 1)) * pip ** (1 / k)
    sqrt_ix = (k + 1) / (k - 1) * (1 - pip ** ((k - 1) / k))
    At_Ae = first_ix * (sqrt_ix) ** 0.5
    return 1 / At_Ae

def th_area(self):
    Pc = self.ch_pressure * 100000; y = self.gamma
    Tc = self.ch_temp; mdot = self.mdot;
    R = 8.314/self.molar_mass
    first = mdot / Pc; second = (R * Tc) / y; third = (2 / (y + 1)); fourth = (y + 1) / (y - 1)
    return first * (second / (third ** fourth)) ** 0.5