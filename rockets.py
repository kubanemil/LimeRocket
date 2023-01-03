from engine import Engine
from tanks import Tanks
from flight import Flight
from tools import *


class Rocket(Flight):
    def __init__(self, tanks_height=1, radius=0.07, ho2_per=38, others_mass=4, mdot=1, o_f=2, ch_pressure=3.5, Po_Pf=2.5, Cd=0.7):
        self.gap_width = 5 * (10**(-3))
        # ratio of oxidizer to fuel, mass flow rate of the propellant, total height of ox+fuel tanks.
        self.o_f, self.mdot, self.radius, self.Cd, self.tanks_height, self.others = o_f, mdot, radius, Cd, tanks_height, others_mass
        self.ho2_dens = 1_450_000/((ho2_per*-4.50)+1450)
        self.ho2_per = ho2_per  # perentage of h2o2 in oxidizer
        self.engine = Engine(ch_pressure=ch_pressure, mdot=mdot, radius=radius-self.gap_width, o_f=self.o_f, ho2_per=ho2_per)     # engine of the rocket
        self.tanks = Tanks(radius=radius-self.gap_width, tanks_height=tanks_height, o_f=self.o_f,
                           mdot=mdot, Po_Pf=Po_Pf, ox_dens=self.ho2_dens)  # oxidizer, gas, fuel tanks
        self.dens = 7900  # kg/m3  # density of the construction steel
        self.burntime = self.tanks.prop_mass / self.mdot  # sec
        self.k = (self.Cd * 1.2 * pi / 2) * (self.radius**2)   # drag coefficient
        self.thick = 0.001  # meters  # thickness of the wall of the tubes.
        self.cap_h = 0.2  # height of the rocket cap
        self.max_mass = 50  # kg maximum allowed mass for rocket
        self.fuelOx_h = 0.05  # distance between fuel and oxidizer tanks
        self.engineGas_h = 0.05  # distance between engine and gas tank
        self.OxCap_h = 0.05  # dist. between oxidizer tank and cap of the rocket
        self.gasFuel_h = 0.15  # dist. b. gas tank and fuel tank
        self.base_h = (self.OxCap_h + self.tanks.gas_height + self.gasFuel_h + self.fuelOx_h + self.tanks.total_height
                       + self.engineGas_h + self.engine.total_length) # height of the base tube
        self.ox_tube = self.tanks.ox_height + (self.OxCap_h+self.fuelOx_h)/2  # height of the oxidizer tube
        self.fuel_tube = self.tanks.fuel_height + (self.gasFuel_h + self.fuelOx_h)/2  # height of the fuel tube
        self.gas_tube = self.tanks.gas_height + (self.engineGas_h + self.gasFuel_h)/2  # height of the gas tube
        self.engine_tube = self.engine.total_length + (self.engineGas_h/2)  # height of the engine tube
        self.fings_h = self.base_h/4    # height of the fings
        self.right_h = self.fings_h/2   # height of the right part
        self.fings_w = self.radius*3  # width of the fings
        self.fings_num = 4  # number of fings
        self.base_mass = (2*pi*self.radius * self.base_h * self.thick)*self.dens
        self.fings_mass = (self.fings_num*self.thick*self.fings_w*(self.fings_h+self.right_h)/2)*self.dens
        self.cap_mass = ((pi*self.radius * ((self.cap_h**2 + self.radius**2)**0.5))*self.thick)*self.dens
        self.total_mass = self.base_mass + self.cap_mass + self.fings_mass + self.engine.total_mass() + \
                            self.tanks.total_mass + self.others
        self.check()

    def check(self):
        if self.total_mass > self.max_mass:
            print(round(self.total_mass, 2), "kg. ", f"Rocket heavier than {self.max_mass} kg")
            return False
        if self.total_mass >= self.engine.thrust/9.8:
            # raise ValueError("Engine can't lift the rocket!")
            print(round(self.total_mass, 2), "kg.", round(self.engine.thrust/9.8, 2), "kg. ", "Engine can't lift it!")
            print(round(self.radius*100, 2), "cm.")
            return False
        return True

    def conf(self):
        """
        Returns the configurations for plotting the rocket.
        """
        return {
            "length": {
                "base": scaled_list([0, 0, self.base_h, self.cap_h]),
                "fings": [0, 0, self.right_h, self.fings_h],
                "engine": self.engine.reversed_lengths,
                "gas_tank": scaled_list([self.engine.total_length + self.engineGas_h,
                                         0, self.tanks.gas_height, 0]),
                "fuel_tank": scaled_list([self.engine.total_length+self.engineGas_h+self.tanks.gas_height+self.gasFuel_h,
                                          0, self.tanks.fuel_height, 0]),
                "ox_tank": scaled_list([self.engine.total_length+self.engineGas_h + self.tanks.gas_height+self.gasFuel_h+
                                         self.tanks.fuel_height + self.fuelOx_h, 0, self.tanks.ox_height, 0]),
                "tubes": scaled_list([0, 0, self.engine_tube, 0, 0, self.gas_tube, 0, 0, self.fuel_tube, 0, 0, self.ox_tube, 0, 0])
            },

            "radius": {
                "base": [0, self.radius, self.radius, 0],
                "fings": [0, self.fings_w + self.radius, self.fings_w+self.radius, self.radius],
                "engine": self.engine.reversed_radiuses,
                "tank": [0, self.tanks.radius, self.tanks.radius, 0],
                "tubes": [0, self.radius, self.radius, 0, self.radius, self.radius, 0, self.radius,
                          self.radius, 0, self.radius, self.radius, 0, self.radius]
            }
        }

    def plot_rocket(self, ix=999, save=False):
        """
        Plots the rocket.
        """
        l_r = self.conf()
        rads = l_r['radius']
        lens = l_r['length']

        plt.rcParams["figure.figsize"] = [13.9, 6.5]
        plt.rcParams["figure.autolayout"] = True
        plt.figure(ix)
        plt.axis([0, 280, -40, 40])
        plt.axis("scaled")
        plot_part(lens['fings'], rads['fings'], fill_color="white")
        plot_part(lens['base'], rads['base'], fill_color="white")
        plot_part(lens['engine'], rads['engine'], fill_color="white")
        plot_part(lens['fuel_tank'], rads['tank'], fill_color="white")
        plot_part(lens['ox_tank'], rads['tank'], fill_color="white")
        plot_part(lens['gas_tank'], rads['tank'], fill_color="white")
        plot_part(lens['tubes'], rads['tubes'], fill_color="white")
        plt.grid(True)
        plt.xlabel("Length, cm")
        plt.ylabel("Radius, cm")
        plt.figtext(0.06, 0.89, s=f"Thrust={round(self.engine.thrust/9.8, 2)} kg", fontsize='12')
        plt.figtext(0.06, 0.85, s=f"ISP={round(self.engine.isp, 2)} sec", fontsize='12')
        plt.figtext(0.21, 0.89, s=f"Total mass={round(self.total_mass, 2)} kg", fontsize='12')
        plt.figtext(0.21, 0.85, s=f"Propellant mass={round(self.tanks.prop_mass, 2)} kg", fontsize='12')
        plt.figtext(0.38, 0.89, s=f"Tanks' height={round(self.tanks_height, 2)} m", fontsize='12')
        plt.figtext(0.38, 0.85, s=f"Burn time={round(self.burntime, 2)} sec", fontsize='12')
        plt.figtext(0.55, 0.89, s=f"Po_Pf={round(self.tanks.Po_Pf, 2)}", fontsize='12')
        plt.figtext(0.65, 0.65, s=self.__repr__(), fontsize='10')
        plt.figtext(0.85, 0.7, s=self.inputs(), fontsize='10')
        dx = 0.2
        gas_x = (lens["gas_tank"][0] + dx)/2.8
        ox_x = (lens["ox_tank"][0] + dx)/2.8
        fuel_x = (lens["fuel_tank"][0] + dx)/2.8
        plt.figtext(gas_x, 0.49, s=f"Gas", fontsize='10')
        plt.figtext(ox_x, 0.49, s=f"OX", fontsize='10')
        plt.figtext(fuel_x, 0.49, s=f"Fuel", fontsize='10')
        plt.figtext(0.07, 0.49, s=f"Engine", fontsize='9')

        if save:
            plt.savefig(f"rocket_plots\\rocket{ix}.png")
        else:
            plt.show()

    def steel_mass(self):
        return self.tanks.tanks_mass() + self.engine.total_mass() + self.cap_mass + self.base_mass + self.fings_mass + 1

    def __repr__(self):
        return f"Rocket thrust: {rnd(self.engine.thrust)} Newtons" "\n" \
               f"ISP: {rnd(self.engine.isp)} secs" "\n" \
               f"Total mass: {rnd(self.total_mass)} kg" "\n" \
               f"Propellant mass: {rnd(self.tanks.prop_mass)} kg" "\n" \
               f"Tank's height: {rnd(self.tanks.total_height)} m" "\n" \
               f"Burn time: {rnd(self.burntime)} secs" "\n" \
               f"Gas pressure {rnd((self.engine.ch_pressure * self.tanks.Po_Pf) + (self.tanks.p_diff_fuel/100000))} bar" "\n" \
               f"Mass flow rate: {rnd(self.mdot)} kg/s" "\n" \
               f"Pressure difference: {rnd(self.tanks.p_diff_fuel/100000)} bar" "\n" \
               f"------------------------" "\n"\
               f"Max Attitude: {rnd(self.max_height())} m" "\n" \

    def inputs(self):
        return f"INPUTS:" "\n" \
            f"Chamber pressure: {rnd(self.engine.ch_pressure)} bar" "\n" \
               f"Radius: {rnd(self.radius)} m" "\n" \
               f"Mdot: {rnd(self.mdot)} kg/s" "\n" \
               f"Other's mass: {rnd(self.others)} kg" "\n" \
               f"Po_Pf: {rnd(self.tanks.Po_Pf)}" "\n" \
               f"Cd: {rnd(self.Cd)}" "\n" \
               f"Tanks height: {rnd(self.tanks_height)} m" "\n" \
               f"O/F ratio: {rnd(self.o_f)}" "\n" \
               f"H2O2 percentage: {rnd(self.ho2_per)}%" "\n"\
               f"H2O2 density: {rnd(self.ho2_dens)} kg/m3" "\n"

if __name__ == "__main__":
    roc = Rocket()
    print(roc)
    # # print("INJ VEL", roc.tanks.fuel_inj_v())
    # # print("OX INJ VEL", roc.tanks.ox_inj_v())
    # # print()
    # roc.plot_rocket()
    print(roc.tanks.p_diff_fuel, roc.tanks.p_diff_ox)

