def At():
    Pc = 2039000
    y = 1.2
    Tc = 2800
    M = 0.02309
    mdot = 1.754
    R = 8.314/M
    first = mdot / Pc
    second = (R*Tc)/y
    third = (2/(y+1))
    fourth = (y+1)/(y-1)
    return first * (second/(third**fourth))**0.5

print(At()*10000)