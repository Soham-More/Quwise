import numpy as np
import matplotlib.pyplot as plt

def jacobi_smoothen(x, y):
    z = np.zeros(len(x))
    for i in range(len(x)):
        if i > 0:
            z[i] += x[i - 1] * 0.5
        if i < len(y) - 1:
            z[i] += x[i + 1] * 0.5
        z[i] += y[i] * 0.5
    return x

def evaluatePoission(x):
    y = np.zeros(len(x))
    for i in range(len(x)):
        if i > 0:
            y[i] += x[i - 1] * -1.0
        if i < len(x) - 1:
            y[i] += x[i + 1] * -1.0
        y[i] += x[i] * 2.0
    return y

def restrict(x, y):
    pass

def multigrid_v_cycle(x, y):
    pass



