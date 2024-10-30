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
    x_unscaled = np.zeros(len(y), dtype=np.float64)

    for i in range(len(y)):
        c_i = 1 - (1 / (i + 2))
        weighted_sum += (i + 1) * y[i]
        x_unscaled[i] = weighted_sum * c_i / ((i + 1)*(i + 1))

    psum = 0
    x = np.zeros(len(y), dtype=np.float64)

    for i in range(len(y)-1, -1, -1):
        psum += x_unscaled[i]
        x[i] = psum * (i + 1)

    return x

def double_direct1(y):
    x = np.array(solvePossion(y))
    if len(y) <= 8:
        return x
    
    r = evaluatePoission(x) - y

    x_r = double_direct1(r[::2])
    x_r = np.interp( np.arange(0, len(y), 1), np.arange(0, len(y), 2), x_r )

    x = x - x_r

    return jacobi_smoothen(x, y)

def double_direct(y):
    x = solvePossion(y)
    r = evaluatePoission(x) - y
    x_r = solvePossion(r)
    x = x - x_r
    return x

def get_eigenvector(i, n):
    return 2 * np.sqrt(2 / (n + 1)) * np.sin( (i * np.arange(1, n + 1, 1, dtype=np.float128) * np.pi) / (n + 1) , dtype=np.float256)

def get_eigenvalue(i, n):
    return 2 * (1 - np.cos( i * np.pi / (n + 1), dtype=np.float256))

N = 2 << 16
ei = N - 1

mask = np.abs(np.linspace(-1.0, 1.0, N)) < 0.1

test = np.zeros(N)
test[mask] = 1.0

test = np.zeros(N) #test / ((N + 1) * (N + 1)) # get_eigenvector(ei, N) # np.random.random(N).astype(np.float128) / (N * N) #np.zeros((N), dtype=np.float128)
test[N//2] = 1.0

test[N//4] = 1.0

test[N//4+N//2] = 1.0

x = solvePossion(test).astype(dtype=np.float64)

y = np.zeros(N, dtype=np.float64)

for i in range(N):
    if i > 0:
        y[i] += x[i - 1] * -1.0
    if i < N - 1:
        y[i] += x[i + 1] * -1.0
    y[i] += x[i] * 2.0

print(np.max(np.abs(y - test)))
print(np.average(np.abs(y - test)))

#print(100 * np.linalg.norm( get_eigenvalue(ei, N) * x - test ) / np.linalg.norm(get_eigenvalue(ei, N) * x))

plt.plot(np.arange(1, N + 1, 1), x / (N + 1), label='x')
plt.plot(np.arange(1, N + 1, 1), test, label='ref')
plt.legend()
plt.show()
