#this unit tests to see if the energy is conserved

import pytest
import numpy as np
from leapfrog import leapfrog


omega1_init= 2

omega2_init=6

theta1_init=np.pi*.8

theta2_init=np.pi

t_final=8

tSteps=300

g=9.8 #m/s^2
L=1 #m
m=1#kg

def test_leapfrog(): 
    omega1List, omega2List, theta1List, theta2List, tList=leapfrog(omega1_init, omega2_init, theta1_init, theta2_init, t_final, tSteps)


    x1=np.sin(theta1List)
    y1=np.cos(theta1List)

    x2=x1+np.sin(theta2List)
    y2=y1+np.cos(theta2List)

    L1=x1**2+y1**2

    L2=(y2-y1)**2+(x2-x1)**2

    print('At each time, step, the legnth of the second pendulum is:')
    print(L2)

    print('number of time steps:')

    print(len(L1))

    assert(L1.all==np.ones(len(L1)).all)
    assert(L2.all==np.ones(len(L2)).all)
test_leapfrog()