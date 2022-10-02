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


    def __repr__(self):
        return f"The maximum height with drag: {self.drag_max_height()/1000} km.\n" + \
            f"The maximum height in drag-free environment: {self.max_height()} km."

# flight = Flight(rocket=Rocket(tanks_height=1))
tank_radiuses = [i for i in range(1, 15)]
drag_flights = [Flight(Rocket(radius=r/100, mdot=1)).drag_max_height() for r in tank_radiuses]
dragfree_flights = [Flight(Rocket(radius=r/100, mdot=1)).max_height() for r in tank_radiuses]
plt.plot(tank_radiuses, drag_flights, label='with drag force')
plt.plot(tank_radiuses, dragfree_flights, label='drag free')
plt.legend()
plt.xlabel("Rocket radius, cm")
plt.ylabel("Max height, meters")

plt.show()

