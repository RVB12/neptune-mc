import numpy as np
import matplotlib.pyplot as plt

def time_since_last_event_exact(T, lamb):
    t = 0
    n = 0
    while True:
        t_new = t + np.random.exponential(1/lamb)
        if t_new > T:
            return T - t, n
        n+=1
        t = t_new

accept_rate = 0
def time_since_last_event_fast(T, lamb):
    # Step 1, sample the number of events
    n = np.random.poisson(lamb * T)
    if(n == 0):
        return T, 0

    # Step 2, sample the time after n events using truncated erlang distribution
    # (truncated because we know t < T)
    # this is sampled using accept-reject
    while True:
        global accept_rate
        accept_rate += 1
        t = np.random.gamma(n, 1/lamb)
        # If t > T, we reject, this is obvious
        if t > T:
            continue
        # We also sample an exponential random variable with rate lambda
        # and check if t + t2 > T, this essentially verifies whether
        # the event is truly the last event
        t2 = np.random.exponential(1/lamb)
        if t + t2 > T:
            return T - t, n



n_samples = 100000

T = 1
lamb = 10

times = np.zeros(n_samples)
n_events = np.zeros(n_samples)

for i in range(n_samples):
    time, n = time_since_last_event_exact(T, lamb)
    times[i] = time
    n_events[i] = n

times2 = np.zeros(n_samples)
n_events2 = np.zeros(n_samples)
for i in range(n_samples):
    time, n = time_since_last_event_fast(T, lamb)
    times2[i] = time
    n_events2[i] = n

print("Accept rate: ", n_samples / accept_rate)

# # Count the contribution for 2 events
# chance = lamb**2 * np.exp(-lamb * T) / 2
# prob += lamb**2 *  np.exp(-lamb * t_plot) * chance

plt.hist(times, bins=30, density=True)
plt.ylabel('Density')
plt.xlabel('Final free flight time')
plt.title(r'$\epsilon = 0.1$')
plt.show()

# fig, ax = plt.subplots(2, 2)
# ax[0,0].hist(times, bins=30, density=True)
# ax[0,0].set_ylabel('Density')
# ax[0,0].set_xlabel('Time since last event (direct)')
# ax[0,1].hist(n_events, bins=30, density=True)
# ax[0,1].set_ylabel('Density')
# ax[0,1].set_xlabel('Number of events')

# ax[1,0].hist(times2, bins=30, density=True)
# ax[1,0].set_ylabel('Density')
# ax[1,0].set_xlabel('Time since last event (fancy sample)')
# ax[1,1].hist(n_events2, bins=30, density=True)
# ax[1,1].set_ylabel('Density')
# ax[1,1].set_xlabel('Number of events')
# plt.show()

# fig, ax = plt.subplots(1, 1)

# counts, bin_edges = np.histogram(times, bins=30, density=True)
# bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
# ax.plot(bin_centers, counts, label='Direct')

# counts, bin_edges = np.histogram(times2, bins=30, density=True)
# bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
# ax.plot(bin_centers, counts, label='Fancy sample')

# plt.show()