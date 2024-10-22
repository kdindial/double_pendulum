# generic rk4 algorithm:
def rk4_general(f, y0, t0, dt): 
  f1 = f(y0, t0)
  f2 = f(y0 + f1*dt/2, t0 + dt/2)
  f3 = f(y0 + f2*dt/2, t0 + dt/2)
  f4 = f(y0 + f3*dt, t0 + dt)
  xout = y0 + (f1 + 2*f2 + 2*f3 + f4)*dt/6
  return xout


def rk4_step_omega1(f, y, t, dt, omega1, omega2, theta1, theta2):
    # Use Runge-Kutta method to update each step for omega1 (f1)
    k1 = f(omega1, omega2, theta1, theta2)
    k2 = f(omega1 + 0.5 * k1 * dt, omega2, theta1, theta2)
    k3 = f(omega1 + 0.5 * k2 * dt, omega2, theta1, theta2)
    k4 = f(omega1 + k3 * dt, omega2, theta1, theta2)

    return y + (k1 + 2*k2 + 2*k3 + k4) * dt / 6


def rk4_step_omega2(f, y, t, dt, omega1, omega2, theta1, theta2):
    # Use Runge-Kutta method to update each step for omega2 (f2)
    k1 = f(omega1, omega2, theta1, theta2)
    k2 = f(omega1, omega2 + 0.5 * k1 * dt, theta1, theta2)
    k3 = f(omega1, omega2 + 0.5 * k2 * dt, theta1, theta2)
    k4 = f(omega1, omega2 + k3 * dt, theta1, theta2)

    return y + (k1 + 2*k2 + 2*k3 + k4) * dt / 6


def rk4(omega1_init, omega2_init, theta1_init, theta2_init, t_final, tSteps):
    #this is a rather faulty implementation of Rungekutta to solve this problem. Id like to keep it here for the time being as we work through this project. 
    # For a better impletmentaion of RK4 see the funciton below labeled "rk4_careful()" 
  # a limitation of this function is that theta1 and theta2 are iterated using euler steps which lack the precision of RK4
  
    t0 = 0
    deltaT = t_final / tSteps

    # time array
    tList = np.arange(t0, t_final, deltaT)

    # store results
    omega1List = []
    omega2List = []
    theta1List = []
    theta2List = []

    # Initialize
    omega1 = omega1_init
    omega2 = omega2_init
    theta1 = theta1_init
    theta2 = theta2_init

    
    for t in tList:
        # add updates to lists
        omega1List.append(omega1)
        omega2List.append(omega2)
        theta1List.append(theta1)
        theta2List.append(theta2)

        # use rk4_step update 
        omega1 = rk4_step_omega1(f1, omega1, t, deltaT, omega1, omega2, theta1, theta2)
        omega2 = rk4_step_omega2(f2, omega2, t, deltaT, omega1, omega2, theta1, theta2)
        theta1 += omega1 * deltaT  # Update angles directly using angular velocity
        theta2 += omega2 * deltaT

    return omega1List, omega2List, theta1List, theta2List, tList



def rk4_careful(omega1_init, omega2_init, theta1_init, theta2_init, t_final, tSteps):
  #this function is better than the previous one because there is a unique k for each variable: omega1, omega2, theta1, theta2

  t0 = 0

    deltaT = t_final/tSteps
    
    tList = np.arange(t0, t_final, step=deltaT)
    omega1List = []
    omega2List = []
    theta1List = []
    theta2List = []

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

        #RK4 updates notice
        #K1
        k1_omega1 = deltaT * f1(omega1, omega2, theta1, theta2)
        k1_omega2 = deltaT * f2(omega1, omega2, theta1, theta2)
        k1_theta1 = deltaT * omega1
        k1_theta2 = deltaT * omega2
        #k2
        k2_omega1 = deltaT * f1(omega1 + k1_omega1/2, omega2 + k1_omega2/2, theta1 + k1_theta1/2, theta2 + k1_theta2/2)
        k2_omega2 = deltaT * f2(omega1 + k1_omega1/2, omega2 + k1_omega2/2, theta1 + k1_theta1/2, theta2 + k1_theta2/2)
        k2_theta1 = deltaT * (omega1 + k1_omega1/2)
        k2_theta2 = deltaT * (omega2 + k1_omega2/2)
        #K3
        k3_omega1 = deltaT * f1(omega1 + k2_omega1/2, omega2 + k2_omega2/2, theta1 + k2_theta1/2, theta2 + k2_theta2/2)
        k3_omega2 = deltaT * f2(omega1 + k2_omega1/2, omega2 + k2_omega2/2, theta1 + k2_theta1/2, theta2 + k2_theta2/2)
        k3_theta1 = deltaT * (omega1 + k2_omega1/2)
        k3_theta2 = deltaT * (omega2 + k2_omega2/2)
        #K4
        k4_omega1 = deltaT * f1(omega1 + k3_omega1, omega2 + k3_omega2, theta1 + k3_theta1, theta2 + k3_theta2)
        k4_omega2 = deltaT * f2(omega1 + k3_omega1, omega2 + k3_omega2, theta1 + k3_theta1, theta2 + k3_theta2)
        k4_theta1 = deltaT * (omega1 + k3_omega1)
        k4_theta2 = deltaT * (omega2 + k3_omega2)

        #new values are Initial value + weighted average of the updates 
        omega1 = omega1 + (k1_omega1 + 2*k2_omega1 + 2*k3_omega1 + k4_omega1)/6
        omega2 = omega2 + (k1_omega2 + 2*k2_omega2 + 2*k3_omega2 + k4_omega2)/6
        theta1 = theta1 + (k1_theta1 + 2*k2_theta1 + 2*k3_theta1 + k4_theta1)/6
        theta2 = theta2 + (k1_theta2 + 2*k2_theta2 + 2*k3_theta2 + k4_theta2)/6

    return omega1List, omega2List, theta1List, theta2List, tList

