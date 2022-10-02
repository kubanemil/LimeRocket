import numpy as np
from rockets import Rocket
import matplotlib.pyplot as plt

class Flight:
    def __init__(self, rocket):
        self.rocket = rocket
        self.engine = rocket.Engine()


    def max_height(self):
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
        return f"The maximum height: {self.max_height()/1000} km."

flight = Flight(rocket=Rocket(fuel_height=0.1))
fuel_heights = [i/10 for i in range(10)]
flights = [Flight(Rocket(fuel_height=h, radius=0.3)).max_height() for h in fuel_heights]
plt.plot(fuel_heights, flights)
plt.show()
print(flight)
