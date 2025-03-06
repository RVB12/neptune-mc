import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-1, 1, 400)
x_mean = 0
sigma_x = 0.1

v = np.linspace(-1, 1, 400)
v_mean = 0
sigma_v = 0.1

X, V = np.meshgrid(x, v)

f0 = np.exp(-((X - x_mean)**2 / (2 * sigma_x**2) + (V - v_mean)**2 / (2 * sigma_v**2)))

# solution at t = 0
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.plot_surface(X, V, f0)
plt.show()

# solution at t = 3
t = 2
f1 = np.exp(-((X - V*t - x_mean - 0.1)**2 / (2 * sigma_x**2) + (V - v_mean)**2 / (2 * sigma_v**2)))

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.plot_surface(X, V, f1)
plt.show()