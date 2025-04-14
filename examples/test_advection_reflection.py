import numpy as np
import matplotlib.pyplot as plt

# Parameters
mu = 5.0      # drift
x0 = 4.4      # initial position
L = 5.0       # reflecting boundary at x = L
t = 0.2       # total time

n = 201
x_grid = np.linspace(0, L, n)
dx = x_grid[1] - x_grid[0]

dt = 0.001
num_steps = int(t / dt)
r2 = mu * dt / dx

density_ghosts = np.zeros(n + 2)
x_index = int(x0 / L * (n)) + 1
density_ghosts[x_index] = 1/dx
# matrix A for implicit method
# assumes positive drift for upwind discretization
# and reflecting boundary at x = L
A = np.zeros((n+2, n+2))
A[0, 1] = 1.0 # neumann condition, ghost cell is equal to the second cell (second order derivative)
A[-1, -1] = 1.0 # neumann condition, ghost cell is equal to the second cell (second order derivative)
for i in range(1, n+1):
    if(i > 1):
        A[i, i] = 1.5 *r2
        A[i, i-1] = - 2 * r2
        A[i, i-2] = 0.5 * r2
    else:
        A[i, i] = r2
        A[i, i-1] = - r2

# now actually solve the system
A1 = np.eye(n+2) - 0.5 * A
A2 = np.eye(n+2) + 0.5 * A
A1[0, 0] = 0.5
A1[-1, -1] = 0.5
A2[0, 0] = -0.5
A2[-1, -1] = 0.5

print(A1)

for step in range(num_steps):
    density_ghosts = np.linalg.solve(A2, A1 @ density_ghosts)
    if(step % 10 == 0):
        density = density_ghosts[1:-1]
        plt.plot(x_grid, density, label='Finite Volume Solution')
        plt.xlabel('x')
        plt.ylabel('Probability Density Function')
        plt.title('Finite volume solution, t = %.3f' % (step * dt))
        plt.legend()
        plt.grid()
        plt.show()

