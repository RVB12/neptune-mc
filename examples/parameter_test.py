import numpy as np
import matplotlib.pyplot as plt

sigma = 1.0
T = 1.0
theta = np.logspace(-4, 4, 100)

D_kinetic = T * np.sqrt(theta) / 2

D_diffusive = T / sigma

D_kd = np.exp(-theta * sigma / 10) * D_kinetic + (1 - np.exp(-theta * sigma / 10)) * D_diffusive


plt.figure()
plt.loglog(theta * sigma, np.sqrt(2 * D_diffusive * theta), label="diffusive limit")
plt.loglog(theta * sigma, np.sqrt(2 * D_kinetic * theta), label="kinetic limit")
plt.loglog(theta * sigma, np.sqrt(2 * D_kd * theta), label="kd term")
plt.xlabel(r"$\theta \sigma$")
plt.legend()
plt.show()