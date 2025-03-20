import numpy as np
import matplotlib.pyplot as plt
import time

def time_since_last_event_exact(T, lamb):
    t = 0
    n = 0
    while True:
        t_new = t + np.random.exponential(1/lamb)
        if t_new > T:
            return T - t, n
        n+=1
        t = t_new

def time_since_last_event_fast(T, lamb):
    # Step 1, sample the number of events
    n = np.random.poisson(lamb * T)
    if(n == 0):
        return T, 0
    
    # Step 2, sample the time of the last event
    # from a beta distribution
    u = np.random.beta(n, 1)
    t = u * T
    return T - t, n

n_samples = 100000

T = 1
lamb = 1.5

times = np.zeros(n_samples)
n_events = np.zeros(n_samples)

walltime1 = time.time()
for i in range(n_samples):
    t, n = time_since_last_event_exact(T, lamb)
    times[i] = t
    n_events[i] = n
print("Sampling time exact: ", time.time() - walltime1)


times2 = np.zeros(n_samples)
n_events2 = np.zeros(n_samples)
walltime2 = time.time()
for i in range(n_samples):
    t, n = time_since_last_event_fast(T, lamb)
    times2[i] = t
    n_events2[i] = n

print("Sampling time fast: ", time.time() - walltime2) 

# # Count the contribution for 2 events
# chance = lamb**2 * np.exp(-lamb * T) / 2
# prob += lamb**2 *  np.exp(-lamb * t_plot) * chance

# plt.hist(times, bins=30, density=True)
# plt.ylabel('Density')
# plt.xlabel('Final free flight time')
# plt.title(r'$\epsilon = 0.1$')
# plt.show()

fig, ax = plt.subplots(2, 2)
ax[0,0].hist(times, bins=30, density=True)
ax[0,0].set_ylabel('Density')
ax[0,0].set_xlabel('Time since last event (direct)')
ax[0,1].hist(n_events, bins=30, density=True)
ax[0,1].set_ylabel('Density')
ax[0,1].set_xlabel('Number of events')

ax[1,0].hist(times2, bins=30, density=True)
ax[1,0].set_ylabel('Density')
ax[1,0].set_xlabel('Time since last event (fancy sample)')
ax[1,1].hist(n_events2, bins=30, density=True)
ax[1,1].set_ylabel('Density')
ax[1,1].set_xlabel('Number of events')
plt.show()

# fig, ax = plt.subplots(1, 1)

# counts, bin_edges = np.histogram(times, bins=30, density=True)
# bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
# ax.plot(bin_centers, counts, label='Direct')

# counts, bin_edges = np.histogram(times2, bins=30, density=True)
# bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
# ax.plot(bin_centers, counts, label='Fancy sample')

# plt.show()