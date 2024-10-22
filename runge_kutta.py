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

