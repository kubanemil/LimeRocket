from math import tan, pi

def A_to_R(A):
    return (A/3.14159)**(0.5)

def R_to_A(R):
    return 3.14159 * (R**2)

def A_to_D(A):
    R = (A/3.14159)**(0.5)
    return 2*R

def D_to_A(D):
    return (3.14159 * (D**2))/4


a = tan(15*pi/180)
print(a)