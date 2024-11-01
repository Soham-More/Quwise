import numpy as np
import matplotlib.pyplot as plt

# not neat, i know
import sys
sys.path.insert(0, 'math/')

from pyvisual import PyVi

pyvi = PyVi('v.pyvi')

i = 0

x = [1.3333333333333334e-06,2.6666666666666668e-06,4.3333333333333339e-06,4.6666666666666672e-06,5.0000000000000004e-06,5.3333333333333337e-06,5.6666666666666669e-06,7.3333333333333331e-06]
x = np.array(x)

y = [1.3333333333333336e-06,3.3333333333333325e-07,3.3333333333333325e-07,3.3333333333333325e-07,3.3333333333333325e-07,3.3333333333333325e-07,1.3333333333333336e-06,1.3333333333333336e-06]
y = np.array(y)

dx = np.zeros(len(x))

for i in range(0, len(x) - 1):
    dx[i] = x[i + 1] - x[i]
dx[-1] = 1e-5 - x[-1]

#print(x[1] + y[1], x[2])

print(y)
print(dx)

pyvi.display_all_sections()

plt.show()

