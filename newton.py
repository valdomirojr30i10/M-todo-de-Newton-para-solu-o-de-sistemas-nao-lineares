
import numpy as np
import snl

def G(x):
    y = np.zeros(len(x))
    y[0] = (x[0]**3) + ((x[0]**2)*x[1]) - (x[0]*x[2]) + 6
    y[1] = (np.e**x[0]) + (np.e**x[1]) - x[2]
    y[2] = (x[1]**2) - (2*x[0]*x[2]) - 4
    return y

def J(x):
    M = np.array([[(3*(x[0]**2)) + (2*x[1]*x[0]) - x[2],     x[0]**2,   -x[0], -G(x)[0]],
                  [               np.e**x[0]           ,  np.e**x[1],    -1  , -G(x)[1]],
                  [                  -2*x[2]           ,      2*x[1], -2*x[0], -G(x)[2]]], dtype=np.float64)
    return M 

x0 = np.array([-1,-2,1])


print(snl.newton(J,x0,1e-5,100))