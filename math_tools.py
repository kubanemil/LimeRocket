from math import tan, pi, e, log

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