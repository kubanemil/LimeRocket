import numpy as np
from rockets import Rocket
import matplotlib.pyplot as plt

class Flight:
    def __init__(self, rocket):
        self.rocket = rocket
        self.engine = rocket.Engine()
        self.tank = rocket.Tanks()
        self.burntime = self.tank.prop_mass()/self.rocket.mdot #sec

    def burntime_velocity(self):
        mass_inx = self.rocket.total_mass() / (self.rocket.total_mass() - (self.rocket.mdot * self.burntime))
        return (self.engine.exit_vel * (np.log(mass_inx))) - (9.8 * self.burntime)

    def max_height(self):
        v = self.burntime_velocity()
        ho = v*self.burntime/2
        h1 = (v*v)/(2*9.8)
        return ho+h1

    def drag_max_height(self):
        dmdt = self.rocket.Engine().mdot
        u = self.engine.exit_vel
        k = self.rocket.K()
        M_init = self.rocket.total_mass()
        M_pt = self.rocket.Tanks().total_mass()
        Thrust = dmdt * u
        btime = M_pt / dmdt
        Mb = M_init - (M_pt / 2);
        Mc = M_init - M_pt

        q = ((Thrust - (Mb * 9.8)) / k) ** (1 / 2)
        p = (2 * k * q) / Mb

        e_inx = np.e ** (-(p * btime))
        vt = q * (1 - e_inx) / (1 + e_inx)

        ln_inx = ((Mc * 9.8) + (k * (vt ** 2))) / (Mc * 9.8)
        hc = (Mc / (2 * k)) * np.log(ln_inx)
        return hc

    def plot_tank_height_to_attitude(self, tank_heights=[1, 2, 3],
                                     tank_rads=[i for i in range(1, 15)]):
        mdot = self.rocket.mdot
        rocket_heights = []
        last_rads, masses, prop_masses = [], [], []
        for h in tank_heights:
            rocs = [Rocket(radius=r / 100, mdot=mdot, tanks_height=h) for r in tank_rads if
                    Rocket(radius=r / 100, mdot=mdot, tanks_height=h).check()]
            flights = [Flight(roc).max_height() for roc in rocs]
            rads = [roc.radius for roc in rocs]
            prop_masses.append(round(rocs[-1].Tanks().prop_mass(), 2))
            masses.append(round(rocs[-1].total_mass(), 2))
            rocket_heights.append(round(rocs[-1].base_h, 2))
            last_rads.append(rocs[-1].radius)

            plt.plot(rads, flights, label=f"Tanks height: {h} m.", marker="o")

        info_dict = {"rocket heights": rocket_heights,
                "tank heights": tank_heights,
                "rocket radius": last_rads,
                "masses": masses,
                "propellant masses": prop_masses,
                }
        print(info_dict)

        plt.figtext(0.1, 0.8, f"Mass Flow {mdot}")
        plt.title("Max heights for Drag-Free Flight of the Rocket.")
        plt.legend()
        plt.xlabel("Rocket radius, cm")
        plt.ylabel("Max height, meters")
        plt.show()
        return info_dict

    def __repr__(self):
        return f"The maximum height with drag: {self.drag_max_height()/1000} km.\n" + \
            f"The maximum height in drag-free environment: {self.max_height()} km."


if __name__ == "__main__":
    roc = Rocket(mdot=1, tanks_height=1)
    flight = Flight(roc)
    # flight.plot_tank_height_to_attitude(tank_heights=[0.5, 1, 1.5])
    print(flight.burntime_velocity())