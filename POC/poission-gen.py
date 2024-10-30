import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def jacobi_smoothen(x, y):
    for i in range(len(x)):
        if i > 0:
            x[i] += x[i - 1] * 0.5
        if i < len(x) - 1:
            x[i] += x[i + 1] * 0.5
        x[i] += y[i] * 0.5
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

def solvePossion(y):
    weighted_sum = 0
    x_unscaled = np.zeros(len(y), dtype=np.float32)

    for i in range(len(y)):
        c_i = 1 - (1 / (i + 2))
        weighted_sum += (i + 1) * y[i]
        x_unscaled[i] = weighted_sum * c_i / ((i + 1)*(i + 1))

    psum = 0
    x = np.zeros(len(y), dtype=np.float32)

    for i in range(len(y)-1, -1, -1):
        psum += x_unscaled[i]
        x[i] = psum * (i + 1)

    return x

def solvePossionBVP(y, a, b, h):
    y *= h*h
    y[0] += a
    y[-1] += b
    return solvePossion(y)

def solveStepUniformPoission(y, t_index, h1, h2):
    # assume x[t_index] = y
    t = 0.0
    for i in range(100):
        x1a = solvePossionBVP(np.copy(y[0:t_index]), 0.0, t, h1)
        x1b = solvePossionBVP(np.copy(y[t_index + 1:]), t, 0.0, h2)
        #t = y[t_index]*h1*h2 + (h2*x1a[-1] + h1*x1b[0]) / (h1 + h2)

        r1 = h2 / (h1 + h2)
        r2 = h1 / (h1 + h2)

        t = ( y[t_index]*h1*h2 + r1 * y[t_index - 1] * h1 * h1 + r2 * y[t_index + 1] * h2 * h2 + r1 * x1a[-2] + r2 * x1b[1] )


    return np.concatenate((x1a, np.array([t]), x1b))

a = 0.0
b = -0.605423

test = np.array([-81.244753, -121.377436, -141.115496, -150.822979, -155.597248, -157.945291, -159.100087, -159.668029, -159.947352, -160.084731, -160.152306, -160.185561, -160.201958, -160.210108, -160.214292, -160.216704, -160.218613, 160.214303, 160.210131, 160.202004, 160.185654, 160.152495, 160.085115, 159.948132, 159.669615, 159.103311, 157.951848, 155.610577, 150.850071, 141.170541, 121.489192, 81.471296])
# np.ones(N)  #np.random.random(N).astype(np.float128) #np.zeros((N), dtype=np.float128)

#test = np.concatenate([np.array([0]*1024), test, np.array([0]*1024)])

N = len(test)

coords = np.linspace(0.0, 2e-5, N + 2)

#x = solvePossionBVP(np.copy(test), a, b, 1/N)
x = solvePossionBVP(np.copy(test), a, b, 1.0/(N+2))

A = np.zeros((N, N))

for i in range(N):
    A[i, i] = 2.0
    if i > 0:
        A[i, i - 1] = -1.0
    if i < N - 1:
        A[i, i + 1] = -1.0

#x = np.linalg.solve(A, test)

x = np.insert(x, 0, a)
x = np.append(x, b)

test = np.insert(test, 0, 0.0)
test = np.append(test, 0.0)

plt.plot(coords, x, label='potential')
#plt.plot(coords, test, label='charge')
plt.legend()
plt.show()
