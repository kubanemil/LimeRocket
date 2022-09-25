""""This graph illustristrates the height to which
it will be elevated. """
import matplotlib.pyplot as plt
import  numpy as np

def K(R, dK=0.75):
    return 0.5*1.2*3.14*(R**2) * dK

# print(K(R=0.06, dK=0.75)*(10**4))

m_release = 12 #kg mass of fuel+LOX
time_release = 12 #sec the time of engine's work
dmdt = m_release/time_release #kg/s
k = 46 *(10**(-4)) #kg/m

u = 2000 #m/s this is specific impulse in space
g = 9.8 #m/s2
Mo = 23 #kg initial mass of the rocket

def terminal_v_rocket(t):
    p1 = (dmdt/k) * (u+(g*t))
    p2 = (Mo*g/k)
    # print(p1,p2)
    return (p1 - p2)**0.5

def time_release_to_reach_terminalV():
        return (Mo/dmdt) - (u/g)

tr = time_release_to_reach_terminalV()
# print(tr) #Time of engine's work, to reach the terminal velocity for give engine's force
# for t in range(0, 100):
#     print(terminal_v_rocket(t))


def height_free_fall(m, k, t):
    cos_inx = t/((m/(g*k))**0.5)
    return (m/k) * np.log(np.cosh(cos_inx))

# print(height_free_fall(m=0.007, k=0.0009, t=0.93))