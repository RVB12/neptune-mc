import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance

def sample_inverse_gaussian(mean, lam):
    """Sample from an inverse Gaussian distribution with given mean and lambda."""
    v = np.random.normal(0, 1)
    y = v**2
    x = mean + (mean**2 * y) / (2 * lam) - (mean / (2 * lam)) * np.sqrt(4 * mean * lam * y + mean**2 * y**2)
    
    return x if np.random.rand() <= mean / (mean + x) else mean**2 / x

def point_to_line_distance(x0, y0, x1, y1, x2, y2):
    """Compute perpendicular distance from point (x0, y0) to line passing through (x1, y1) and (x2, y2)."""
    return abs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1) / np.sqrt((y2 - y1)**2 + (x2 - x1)**2)

def first_exit_time_inverse_gaussian_triangle(x0, y0, vertices, D):
    """Estimate first exit time from a triangle using the inverse Gaussian method."""
    exit_times = []
    
    for (x1, y1), (x2, y2) in zip(vertices, vertices[1:] + [vertices[0]]):  # Loop over edges
        dist = point_to_line_distance(x0, y0, x1, y1, x2, y2)
        exit_time = sample_inverse_gaussian(dist / D, (dist ** 2) / (2 * D))
        exit_times.append(exit_time)
    
    return min(exit_times)

def is_inside_triangle(x, y, vertices):
    """Check if point (x, y) is inside the triangle using barycentric coordinates."""
    (x1, y1), (x2, y2), (x3, y3) = vertices
    detT = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
    
    alpha = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / detT
    beta  = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / detT
    gamma = 1 - alpha - beta
    
    return (alpha >= 0) and (beta >= 0) and (gamma >= 0)

def first_exit_time_direct_simulation_triangle(x0, y0, vertices, D, dt=0.001):
    """Simulate Brownian motion inside a triangle until it exits."""
    x, y = x0, y0
    t = 0
    sqrt_2Ddt = np.sqrt(2 * D * dt)
    
    while is_inside_triangle(x, y, vertices):
        x += sqrt_2Ddt * np.random.normal()
        y += sqrt_2Ddt * np.random.normal()
        t += dt

    return t

# Define Triangle
triangle_vertices = [(0, 0), (1, 0), (0.0, 0.4)]  # Triangle with vertices (A, B, C)
x0, y0 = 0.5, 0.05  # Start inside the triangle
D = 0.1            # Diffusion coefficient
num_trials = 10000  # Number of simulations

# Run simulations
exit_times_ig = np.array([first_exit_time_inverse_gaussian_triangle(x0, y0, triangle_vertices, D) for _ in range(num_trials)])
exit_times_sim = np.array([first_exit_time_direct_simulation_triangle(x0, y0, triangle_vertices, D) for _ in range(num_trials)])

# Plot the results
plt.hist(exit_times_ig, bins=100, alpha=0.6, color='b', label="Inverse Gaussian Method", density=True, range=(0, 5))
plt.hist(exit_times_sim, bins=100, alpha=0.6, color='r', label="Direct Simulation", density=True, range=(0, 5))
plt.xlabel("Exit Time")
plt.ylabel("Probability Density")
plt.legend()
plt.title("Comparison of Exit Time Distributions for a Triangle")
plt.show()

# Print mean exit times
print(f"Mean Exit Time (Inverse Gaussian): {np.mean(exit_times_ig):.4f}")
print(f"Mean Exit Time (Direct Simulation): {np.mean(exit_times_sim):.4f}")