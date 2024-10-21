import numpy as np

def f1(omega1, omega2, theta1, theta2, g = 9.81, l=1):
    num = (omega1**2)*np.sin(2*theta1 - 2*theta2) + 2*(omega2**2)*np.sin(theta1-theta2) + (g/l)*(np.sin(theta1 - 2*theta2) + 3*np.sin(theta1))
    den = 3 - np.cos(2*theta1 - 2*theta2)

    return -num/den

def f2(omega1, omega2, theta1, theta2, g = 9.81, l=1):
    num = 4*(omega1**2)*np.sin(theta1 - theta2) + (omega2**2)*np.sin(2*theta1-2*theta2) + (2*g/l)*(np.sin(2*theta1 - theta2) - np.sin(theta2))
    den = 3 - np.cos(2*theta1 - 2*theta2)

    return num/den

