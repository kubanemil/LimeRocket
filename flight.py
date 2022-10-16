import numpy as np
from rockets import Rocket
import matplotlib.pyplot as plt

class Flight:
    def __init__(self, rocket):
        self.rocket = rocket
        self.burntime = self.rocket.tanks.prop_mass()/self.rocket.mdot #sec

    def burntime_velocity(self):
        mass_inx = self.rocket.total_mass() / (self.rocket.total_mass() - (self.rocket.mdot * self.burntime))
        return (self.rocket.engine.exit_vel * (np.log(mass_inx))) - (9.8 * self.burntime)

    def dragfree_max_height(self):
        v = self.burntime_velocity()
        ho = v*self.burntime/2
        h1 = (v*v)/(2*9.8)
        return ho+h1

    def max_height(self):
        roc = self.rocket
        k = roc.k()
        m_init = roc.total_mass()
        m_final = roc.total_mass() - roc.tanks.prop_mass()
        T = roc.engine.thrust
        m = m_init
        H, v = 0, 0
        dt = 0.0001
        while True:
            if m > m_final:
                a = self.rocket_acc(T, m, k, v)
                v = v + a*dt
                H += v*dt
                m -= roc.mdot*dt
            if m <= m_final:
                self.max_vel = v
                H += self.elevation_with_drag(v, m_final, k)
                return H

    def plot_tank_height_to_attitude(self, tank_heights=[1, 2, 3],
                                     tank_rads=[i for i in range(1, 15)]):
        mdot = self.rocket.mdot
        rocket_heights = []
        last_rads, masses, prop_masses, max_heights = [], [], [], []
        for h in tank_heights:
            rocs = []
            for r in tank_rads:
                roc = Rocket(radius=r / 100, mdot=mdot, tanks_height=h)
                if roc.check():
                    rocs.append(roc)
                else:
                    break
            # rocs = [Rocket(radius=r / 100, mdot=mdot, tanks_height=h) for r in tank_rads if
            #         Rocket(radius=r / 100, mdot=mdot, tanks_height=h).check()]
            flights = [Flight(roc).max_height() for roc in rocs]
            rads = [roc.radius for roc in rocs]
            prop_masses.append(round(rocs[-1].tanks.prop_mass(), 2))
            masses.append(round(rocs[-1].total_mass(), 2))
            rocket_heights.append(round(rocs[-1].base_h, 2))
            last_rads.append(rocs[-1].radius)
            max_heights.append(round(flights[-1], 2))

            plt.plot(rads, flights, label=f"Tanks height: {h} m.", marker="o")

        info_dict = {
            "attitudes": max_heights,
            "rocket heights": rocket_heights,
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
        plt.close()
        return info_dict

    def __repr__(self):
        return f"The maximum height with drag: {self.max_height()/1000} km.\n" + \
            f"The maximum height in drag-free environment: {self.max_height()} km."

    @staticmethod
    def drag_acc(v, m, k):
        return -(k*(v**2)/m) - 9.8

    @staticmethod
    def rocket_acc(T, m, k, v):
        return T/m - (k/m)*(v**2) - 9.8

    def elevation_with_drag(self, v_init, m, k):
        dt = 0.001
        H, T = 0, 0
        v = v_init
        while True:
            a = self.drag_acc(v, m, k)
            v = v + a*dt
            H += v*dt
            T += dt
            if v <= 0:
                return H


if __name__ == "__main__":
    roc = Rocket(radius=0.06, mdot=0.5, tanks_height=1)
    flight = Flight(roc)
    flight_info = flight.plot_tank_height_to_attitude(tank_heights=[1, 1.5, 2])
    max_height = flight_info["attitudes"][0]


    # print(mdots)
    # print(max_heights)
    # plt.plot(mdots, max_heights)
    # plt.show()