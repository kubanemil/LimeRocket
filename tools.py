from math import tan, pi, e, log
import numpy as np
import matplotlib.pyplot as plt


def octagon_A2R(A):
    return np.sqrt(A/np.sqrt(8))


def octagon_R2A(R):  # radius
    return np.sqrt(8)*R**2

def octagon_R2s(R):
    return R*np.sqrt(2 - np.sqrt(2))
def octagon_s2R(s):
    return s/np.sqrt(2 - np.sqrt(2))
def octagon_A2s(A):
    R = octagon_A2R(A)
    return octagon_R2s(R)

def octagon_s2A(A):
    R = octagon_A2R(A)
    return octagon_R2s(R)

def K2C(K):  # Kelvin to Celsius
    return K - 273

def to_mm(m):
    return m*1000

def rnd(f):
    return round(f, 3)


def A2R(A):
    return (A/3.14159)**(0.5)


def R2A(R):
    return 3.14159 * (R**2)


def A2D(A):
    R = (A/3.14159)**(0.5)
    return 2*R


def D2A(D):
    return (3.14159 * (D**2))/4


def scaled_list(lst):
    new_lst = []
    for i in range(len(lst)):
        if i == 0:
            new_lst.append(lst[i])
        else:
            new_lst.append(lst[i] + new_lst[i - 1])
    return new_lst


def to_bar(Pascal):
    return Pascal/(10**5)


def plot_part(length_list, radius_list, fill_color="white", line_color="black"):
    lengths_cm = np.array(length_list) * 100
    radiuses_cm = np.array(radius_list) * 100
    plt.plot(lengths_cm, radiuses_cm, line_color, linewidth=0.5)
    plt.plot(lengths_cm, -radiuses_cm, line_color, linewidth=0.5)
    plt.fill_between(lengths_cm, radiuses_cm, color=fill_color)
    plt.fill_between(lengths_cm, -radiuses_cm, color=fill_color)

