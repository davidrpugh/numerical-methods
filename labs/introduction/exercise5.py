from __future__ import division
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

# set params
N = 100
T = 10000
p = 0.75
loaded_coin = stats.distributions.bernoulli(p)

# generate an array of integers
data = loaded_coin.rvs((T,N))

# create an array in which to store our sample averages
loaded_difference = np.zeros(data.shape, dtype=float)

for j in range(N):
    for i in range(T):
        loaded_difference[i, j] = 2 * np.sum(data[:i + 1, j]) - (i + 1)
        
# create new Figure and Axes objects
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)

# plot each sample path...
for i in range(N):
    ax.plot(np.arange(1, T + 1, 1), loaded_difference[:, i], 'k-', alpha=0.05)

# label axes
xlab = ax.set_xlabel('Index', fontsize=15)
xlab.set_family('serif')
ylab = ax.set_ylabel('Heads - Tails', family='serif', fontsize=15)

# set the title
title = ax.set_title('Difference really diverges with a loaded coin!')
title.set_family('serif')
title.set_fontsize(15)

plt.savefig('graphics/exercise-5a.png')
plt.show()

# set params
N = 100
T = 10000

# generate an array of integers 
data = loaded_coin.rvs((T,N))

# create an array in which to store our sample averages
sample_averages = np.zeros(data.shape, dtype=float)

for j in range(N):
    for i in range(T):
        sample_averages[i, j] = np.mean(data[:i + 1, j])
        
# create new Figure and Axes objects
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)

# plot each sample path from obs. 10 forward...
for i in range(N):
    ax.plot(np.arange(1, T + 1, 1), sample_averages[:, i], 'k-', alpha=0.05)
ax.set_ylim(0, 1)

# set the x-axis to have log scale
ax.set_xscale('log')

# demarcate the true mean
ax.axhline(y=p, color='black', linestyle='dashed', label=r'$\mu$')

# label axes
xlab = ax.set_xlabel('Index')
xlab.set_family('serif')

ylab = ax.set_ylabel(r'$\hat{\mu}$')
ylab.set_rotation('horizontal')
ylab.set_fontsize(15)

# set the title
title = ax.set_title(r'Average number of heads converges to $\mu=%g$!' %p)
title.set_family('serif')
title.set_fontsize(15)

# add the legend
ax.legend(loc='best', frameon=False)

plt.savefig('graphics/exercise-5b.png')