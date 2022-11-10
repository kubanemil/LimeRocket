import numpy as np

class Flight:
    def burntime_velocity(self):
        mass_inx = self.total_mass / (self.total_mass - (self.mdot * self.burntime))
        return (self.engine.exit_vel * (np.log(mass_inx))) - (9.8 * self.burntime)

    def dragfree_max_height(self):
        v = self.burntime_velocity()
        ho = v*self.burntime/2
        h1 = (v*v)/(2*9.8)
        return ho+h1

    def max_height(self):
        k = self.k
        m_init = self.total_mass
        m_final = self.total_mass - self.tanks.prop_mass
        T = self.engine.thrust
        m = m_init
        H, v = 0, 0
        dt = 0.0001
        while True:
            if m > m_final:
                a = self.rocket_acc(T, m, k, v)
                v = v + a*dt
                H += v*dt
                m -= self.mdot*dt
            if m <= m_final:
                self.max_vel = v
                H += self.elevation_with_drag(v, m_final, k)
                return H

    def __repr__(self):
        return f"The maximum height with drag: {self.max_height()/1000} km.\n" + \
            f"The maximum height in drag-free environment: {self.dragfree_max_height()/1000} km."

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
