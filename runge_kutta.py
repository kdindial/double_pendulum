#rk4 
def rk4(f, y0, t0, dt): 
  f1 = f(y0, t0)
  f2 = f(y0 + f1*dt/2, t0 + dt/2)
  f3 = f(y0 + f2*dt/2, t0 + dt/2)
  f4 = f(y0 + f3*dt, t0 + dt)
  xout = y0 + (f1 + 2*f2 + 2*f3 + f4)*dt/6
  return xout
