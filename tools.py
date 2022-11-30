from math import tan, pi, e, log
import numpy as np
import matplotlib.pyplot as plt


def to_mm(m):
    return m*1000

def rnd(f):
    return round(f, 3)


def A_to_R(A):
    return (A/3.14159)**(0.5)


def R_to_A(R):
    return 3.14159 * (R**2)


def A_to_D(A):
    R = (A/3.14159)**(0.5)
    return 2*R


def D_to_A(D):
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
