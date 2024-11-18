from omegaDots import f1, f2
import numpy as np
import math

def leapfrog(omega1_init, omega2_init, theta1_init, theta2_init, t_final, tSteps, delta=1e-4):
    # This function implements the leapfrog method of numerical integration
    #
    # INPUTS:
    # omega1_init, omega2_init: initial omegas
    # theta1_init, theta2_init: initial thetas
    # t_final: the final time to which the system is allowed to evolve
    # tSteps: the number of steps taken from the initial time (here assumed to be zero) to f_final
    # delta: the specified error target level
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

    deltaT1Omega = deltaT
    deltaT2Omega = deltaT
    deltaT1Theta = deltaT
    deltaT2Theta = deltaT

    for t in tList:
        # update lists
        omega1List.append(omega1)
        omega2List.append(omega2)
        theta1List.append(theta1)
        theta2List.append(theta2)

        # Now, for each time we evolve using equations 29 of ODE notes
        # first eqn
        k1_1Omega = deltaT1Omega * f1(omega1, omega2, theta1, theta2)
        k1_2Omega = deltaT2Omega * f2(omega1, omega2, theta1, theta2)
        k1_1Theta = deltaT1Theta * omega1
        k1_2Theta = deltaT2Theta * omega2
        # second eqn
        w1HalfOmega = w1HalfOmega + k1_1Omega
        w2HalfOmega = w2HalfOmega + k1_2Omega
        w1HalfTheta = w1HalfTheta + k1_1Theta
        w2HalfTheta = w2HalfTheta + k1_2Theta
        # third eqn
        k2_1Omega = deltaT1Omega * f1(w1HalfOmega, w2HalfOmega, w1HalfTheta, w2HalfTheta)
        k2_2Omega = deltaT2Omega * f2(w1HalfOmega, w2HalfOmega, w1HalfTheta, w2HalfTheta)
        k2_1Theta = deltaT1Theta * w1HalfOmega
        k2_2Theta = deltaT2Theta * w2HalfOmega

        # Now, check the error which will occur, as omega2-omega1 = k2 in section 3 of the ODE notes
        epsilon1Omega = np.abs(k2_1Omega)/30
        epsilon2Omega = np.abs(k2_2Omega)/30
        epsilon1Theta = np.abs(k2_1Theta)/30
        epsilon2Theta = np.abs(k2_2Theta)/30

        if (epsilon1Omega > delta*deltaT1Omega):
            deltaTn1Omega = deltaT1Omega * (30*delta*deltaT1Omega / np.abs(k2_1Omega))**(1/4)
            deltaT1Omega = deltaTn1Omega
            # third eqn
            k2_1Omega = deltaT1Omega * f1(w1HalfOmega, w2HalfOmega, w1HalfTheta, w2HalfTheta)
        
        if (epsilon2Omega > delta*deltaT2Omega):
            deltaTn2Omega = deltaT2Omega * (30*delta*deltaT2Omega / np.abs(k2_2Omega))**(1/4)
            deltaT2Omega = deltaTn2Omega
            # third eqn
            k2_2Omega = deltaT2Omega * f2(w1HalfOmega, w2HalfOmega, w1HalfTheta, w2HalfTheta)

        if (epsilon1Theta > delta*deltaT1Theta):
            deltaTn1Theta = deltaT1Theta * (30*delta*deltaT1Theta / np.abs(k2_1Theta))**(1/4)
            deltaT1Theta = deltaTn1Theta
            # third eqn
            k2_1Theta = deltaT1Theta * w1HalfOmega

        if (epsilon2Theta > delta*deltaT2Theta):
            deltaTn2Theta = deltaT2Theta * (30*delta*deltaT2Theta / np.abs(k2_2Theta))**(1/4)
            deltaT2Theta = deltaTn2Theta
            # third eqn
            k2_2Theta = deltaT2Theta * w2HalfOmega

        # fourth eqn
        omega1 = omega1 + k2_1Omega
        omega2 = omega2 + k2_2Omega
        theta1 = theta1 + k2_1Theta
        theta2 = theta2 + k2_2Theta

        """ 
        I think I need to add the shift to deltaTn here as well for epsilon < junk. Also, need to think about the 2 in t + 2DeltaT in section 3
        """


    return omega1List, omega2List, theta1List, theta2List, tList
