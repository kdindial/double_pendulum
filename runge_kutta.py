import numpy as np
import numpy as np
import math
import matplotlib.pyplot as plt
from omegaDots import f1, f2

def rk4_step(f1, f2, omega1, omega2, theta1, theta2, deltaT):
        """
        Perform a single RK4 step for omega1, omega2, theta1, and theta2.
        """
        # K1
        k1_omega1 = deltaT * f1(omega1, omega2, theta1, theta2)
        k1_omega2 = deltaT * f2(omega1, omega2, theta1, theta2)
        k1_theta1 = deltaT * omega1
        k1_theta2 = deltaT * omega2
        # K2
        k2_omega1 = deltaT * f1(omega1 + k1_omega1 / 2, omega2 + k1_omega2 / 2, theta1 + k1_theta1 / 2, theta2 + k1_theta2 / 2)
        k2_omega2 = deltaT * f2(omega1 + k1_omega1 / 2, omega2 + k1_omega2 / 2, theta1 + k1_theta1 / 2, theta2 + k1_theta2 / 2)
        k2_theta1 = deltaT * (omega1 + k1_omega1 / 2)
        k2_theta2 = deltaT * (omega2 + k1_omega2 / 2)
        # K3
        k3_omega1 = deltaT * f1(omega1 + k2_omega1 / 2, omega2 + k2_omega2 / 2, theta1 + k2_theta1 / 2, theta2 + k2_theta2 / 2)
        k3_omega2 = deltaT * f2(omega1 + k2_omega1 / 2, omega2 + k2_omega2 / 2, theta1 + k2_theta1 / 2, theta2 + k2_theta2 / 2)
        k3_theta1 = deltaT * (omega1 + k2_omega1 / 2)
        k3_theta2 = deltaT * (omega2 + k2_omega2 / 2)
        # K4
        k4_omega1 = deltaT * f1(omega1 + k3_omega1, omega2 + k3_omega2, theta1 + k3_theta1, theta2 + k3_theta2)
        k4_omega2 = deltaT * f2(omega1 + k3_omega1, omega2 + k3_omega2, theta1 + k3_theta1, theta2 + k3_theta2)
        k4_theta1 = deltaT * (omega1 + k3_omega1)
        k4_theta2 = deltaT * (omega2 + k3_omega2)

        # Compute next values
        omega1_next = omega1 + (k1_omega1 + 2 * k2_omega1 + 2 * k3_omega1 + k4_omega1) / 6
        omega2_next = omega2 + (k1_omega2 + 2 * k2_omega2 + 2 * k3_omega2 + k4_omega2) / 6
        theta1_next = theta1 + (k1_theta1 + 2 * k2_theta1 + 2 * k3_theta1 + k4_theta1) / 6
        theta2_next = theta2 + (k1_theta2 + 2 * k2_theta2 + 2 * k3_theta2 + k4_theta2) / 6

        return omega1_next, omega2_next, theta1_next, theta2_next

def rk4_adaptive(omega1_init, omega2_init, theta1_init, theta2_init, t_final, deltaT_init=0.01, tol = 0.00001):

    # Initialize variables
    i = 0
    t = 0
    deltaT = deltaT_init
    omega1, omega2, theta1, theta2 = omega1_init, omega2_init, theta1_init, theta2_init

    # Lists to store results
    tList = [t]
    omega1List = [omega1]
    omega2List = [omega2]
    theta1List = [theta1]
    theta2List = [theta2]


    # main loop 
    while t < t_final:
        i += 1 # index for counting and trouble shooting
      
        # Compute 2 RK4 steps with current deltaT
        mid_omega1, mid_omega2, mid_theta1, mid_theta2 = rk4_step(f1, f2, omega1, omega2, theta1, theta2, deltaT )
        omega1_rk4, omega2_rk4, theta1_rk4, theta2_rk4 = rk4_step(f1, f2, mid_omega1, mid_omega2, mid_theta1, mid_theta2, deltaT )
        # Compute RK4 step with double deltaT
        omega1_rk4_double, omega2_rk4_double, theta1_rk4_double, theta2_rk4_double = rk4_step(f1, f2, omega1, omega2, theta1, theta2, 2 * deltaT)


        # Error estimate (difference between single and double step)
        # because the time step is used for all for updates , only the higheset error value of the 4 variables are used to control the adaptive time step
        error_omega1 = np.abs(omega1_rk4 - omega1_rk4_double)
        error_omega2 = np.abs(omega2_rk4 - omega2_rk4_double)
        error_theta1 = np.abs(theta1_rk4 - theta1_rk4_double)
        error_theta2 = np.abs(theta2_rk4 - theta2_rk4_double)
        max_error = max(error_omega1, error_omega2, error_theta1, error_theta2)
        
        #bad things happen when the errors all = 0 , so when this happens accept the step without calculating rho
        # the deltaT with be changed to 2deltaT
        if max_error == 0:
          rho = 1000
        
        else:
          rho = 30 * deltaT * tol / max_error #eqn 8.53 pg 358

        # Adjust time step based on error
        #if max_error < tol:
        if rho > 1:
            # Accept step
            # advance t by this deltaT *2
            t += 2 * deltaT
            #assign accepted values to lists
            omega1, omega2, theta1, theta2 = omega1_rk4, omega2_rk4, theta1_rk4, theta2_rk4
            tList.append(t)
            omega1List.append(omega1)
            omega2List.append(omega2)
            theta1List.append(theta1)
            theta2List.append(theta2)
            
            #Increase deltaT to  --> h'
      
            deltaT *= min(2,rho ** 0.25) #8.52
            #print(f"{i}up, t={t}, deltaT={deltaT}, max_error={max_error}, rho={rho}")
            
        else:
            # Reject step and deltaT --> h'
            deltaT *= max(.5,rho ** 0.25)
            #print(f"{i}down, t={t}, deltaT={deltaT}, max_error={max_error}, rho={rho}")

    return omega1List, omega2List, theta1List, theta2List, tList,

