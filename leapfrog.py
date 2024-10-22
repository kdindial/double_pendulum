from omegaDots import f1, f2
import numpy as np
import math

# || INITIAL CONDITIONS
theta1_0 = math.radians(0)
theta2_0 = math.radians(3)
omega1_0 = 0 # initially at rest
omega2_0 = 0 # initially at rest

t0 = 0
tf = 100.
steps = 10000

def leapfrog(omega1_init, omega2_init, theta1_init, theta2_init, t_final, tSteps):
    # This function implements the leapfrog method of numerical integration
    #
    # INPUTS:
    # omega1_init, omega2_init: initial omegas
    # theta1_init, theta2_init: initial thetas
    # t_final: the final time to which the system is allowed to evolve
    # tSteps: the number of steps taken from the initial time (here assumed to be zero) to f_final
    #
    # OUTPUTS:
    # omega1List, omega2List: lists of the omegas at each time
    # theta1List, theta2List: lists of the thetas at each time
    # tList: the actual list of times

    t0 = 0

    deltaT = t_final/tSteps
    tList = np.arange(t0, t_final, step=deltaT)

    omega1List = []
    omega2List = []
    theta1List = []
    theta2List = []

    # First step is to do half an Euler step for omega and theta
    w1HalfOmega = omega1_init + (deltaT/2) * f1(omega1_init, omega2_init, theta1_init, theta2_init)
    w2HalfOmega = omega2_init + (deltaT/2) * f2(omega1_init, omega2_init, theta1_init, theta2_init)
    w1HalfTheta = theta1_init + (deltaT/2) * omega1_init
    w2HalfTheta = theta2_init + (deltaT/2) * omega2_init
    

    omega1 = omega1_init
    omega2 = omega2_init
    theta1 = theta1_init
    theta2 = theta2_init

    for t in tList:
        # update lists
        omega1List.append(omega1)
        omega2List.append(omega2)
        theta1List.append(theta1)
        theta2List.append(theta2)

        # Now, for each time we evolve using equations 29 of ODE notes
        # first eqn
        k1_1Omega = deltaT * f1(omega1, omega2, theta1, theta2)
        k1_2Omega = deltaT * f2(omega1, omega2, theta1, theta2)
        k1_1Theta = deltaT * omega1
        k1_2Theta = deltaT * omega2
        # second eqn
        w1HalfOmega = w1HalfOmega + k1_1Omega
        w2HalfOmega = w2HalfOmega + k1_2Omega
        w1HalfTheta = w1HalfTheta + k1_1Theta
        w2HalfTheta = w2HalfTheta + k1_2Theta
        #third eqn
        k2_1Omega = deltaT * f1(w1HalfOmega, w2HalfOmega, w1HalfTheta, w2HalfTheta)
        k2_2Omega = deltaT * f2(w1HalfOmega, w2HalfOmega, w1HalfTheta, w2HalfTheta)
        k2_1Theta = deltaT * w1HalfOmega
        k2_2Theta = deltaT * w2HalfOmega
        # fourth eqn
        omega1 = omega1 + k2_1Omega
        omega2 = omega2 + k2_2Omega
        theta1 = theta1 + k2_1Theta
        theta2 = theta2 + k2_2Theta

    return omega1List, omega2List, theta1List, theta2List, tList