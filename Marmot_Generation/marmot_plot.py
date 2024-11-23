import numpy as np
import matplotlib.pyplot as plt

# Read data
data = np.loadtxt("Marmot/marmot_data.csv", delimiter=",", skiprows=1)

def polynomial(x, a, b, c, d, e):
    return a + b*x + c*x**2 + d*x**3 + e*x**4

def solve(y, a, b, c, d, e):
    return np.roots([e, d, c, b, a-y])

# Know that at x=0, y=3.95 and x=1 y=6.71
# x data is in column 0, y data is in column 1
# x data is not to right scale, and may be shifted
# y data is right scale and not shifted

# Fit a polynomial to the data
print(data[:,0])
from scipy.optimize import curve_fit
popt, pcov = curve_fit(polynomial, data[:,0], data[:,1])

# Plot the data and the fit
plt.plot(data[:,0], data[:,1], 'o', label='data')
x = np.linspace(50,250,100)
plt.plot(x, polynomial(x, *popt), label='fit')
plt.legend()
plt.show()
