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