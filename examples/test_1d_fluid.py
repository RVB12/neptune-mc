import numpy as np
import matplotlib.pyplot as plt

n = 100

t_final = 1.0e-4

x = np.linspace(0, 1, n+1)

x0 = 0.985

ion_density = 1.0e21
ion_velocity = 0
ion_temperature = 1

K_cx = 3.2e-15 * np.sqrt(ion_temperature/0.026)
collision_rate = K_cx * ion_density

BGK_stddev = np.sqrt(ion_temperature * 1.602176634e-19 / 1.67e-27)

print('Collision rate: ', collision_rate)
print('Average collision time: ', 1/collision_rate)

print(x)


D = BGK_stddev**2 / collision_rate
# simulate using forward euler
h = 1 / n
dt = 1.0e-6
r = D * dt / h**2
print('r = ', r)

# this vector includes two ghost cells

density_fluid = np.zeros(n+3)
x_index = int(x0 * (n+1)) + 1
density_fluid[x_index] = n

for t in np.arange(0, t_final, dt):
    # enforce boundary conditions
    density_fluid[0] = density_fluid[1]
    density_fluid[-1] = density_fluid[-2]

    # update density for interior points
    density_fluid[1:-1] = density_fluid[1:-1] + r * (density_fluid[2:] - 2 * density_fluid[1:-1] + density_fluid[:-2])

# remove ghost cells
density_fluid = density_fluid[1:-1]

# simulate using backward euler
h = 1 / n
dt = 1.0e-5
r = D * dt / h**2

density_fluid2 = np.zeros(n+1)
x_index = int(x0 * (n+1))
density_fluid2[x_index] = n

# matrix A for implicit method
A = np.zeros((n+1, n+1)) 
A[0, 0] = 1 + r
A[0, 1] = -r
A[-1, -1] = 1 + r
A[-1, -2] = -r
for i in range(1, n):
    A[i, i] = 1 + 2 * r
    A[i, i-1] = -r
    A[i, i+1] = -r

for t in np.arange(0, t_final, dt):
    density_fluid2 = np.linalg.solve(A, density_fluid2)



# simulate using monte carlo
n_particles = 10000
density_mc = np.zeros_like(x)
for i in range(n_particles):
    x_particle = x0
    t = 0
    v = np.random.normal(ion_velocity, BGK_stddev)

    while t < t_final:
        # time till next collision
        eps = np.random.rand()
        tc = -np.log(eps) / collision_rate
        tc = min(tc, t_final - t)

        # calculate time till boundary hit
        t_boundary = np.inf
        if v > 0:
            t_boundary = (1 - x_particle) / v
        elif v < 0:
            t_boundary = x_particle / -v

        if t_boundary < tc:
            # move particle to boundary
            x_particle += v * t_boundary
            v = -v
            t += t_boundary
            continue

        # move particle
        x_particle += v * tc
        # draw new velocity
        v = np.random.normal(ion_velocity, BGK_stddev)
        t += tc

    # tally final position
    x_index = int(x_particle * (n+1))
    volume = 1/(n)
    density_mc[x_index] += 1/volume / n_particles


# Analytical solution (not valid near boundaries)
density_analytical = np.zeros_like(x)
for i in range(n+1):
    density_analytical[i] = np.exp( -( x[i] - x0 )**2/ (4*D*t_final) ) / np.sqrt(4*np.pi*D*t_final)


plt.plot(x, density_mc, label='Monte Carlo')
plt.plot(x, density_analytical, label='Analytical')
plt.plot(x, density_fluid2, label='Finite Difference')
plt.legend()
plt.xlabel('x')
plt.ylabel('density')
plt.show()