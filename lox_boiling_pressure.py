import matplotlib.pyplot as plt
import numpy as np
from math import e

L = 6820 #J/mol heat of vapor for O2
R = 8.31 # J/K.mol
Pref = 10000 #Pascal
Tref = 90 #Kelvin
def P(T):
    inx = (L/R) * ( (1/Tref)-(1/T) )
    return Pref * (e ** inx)

for i in range(90, 320, 10):
    print(i-273, "Celsius -", P(i)//100000, "Bar")


#can use LOX until it will heat up for +40 from 90 K. Consider using O2 gas.