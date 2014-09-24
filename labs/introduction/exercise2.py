import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

# specify degrees of freedom
df = 3

# create a random variable object representing the normal distribution
chi_square_rv = stats.distributions.chi2(df)

# generate a large sample
T = 10000
sample = chi_square_rv.rvs((T,))

# create new Figure and Axes objects
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)

# histogram the sample
n, bins, patches = ax.hist(sample, color='purple', bins = 25, normed=True, alpha = 0.75)
ax.set_xlim(0, 15)

# overlay the theoretical pdf
grid = np.linspace(0, 15, 1000)
ax.plot(grid, chi_square_rv.pdf(grid), 'r--', label=r'$\chi^2(df)$')

# label the axes
ax.set_xlabel('X', fontsize=15, family='serif')
ax.set_ylabel('Probability density', fontsize=15, family='serif')

# add a title
ax.set_title('We can use the rvs() method to check the pdf() method!', 
             family='serif', fontsize=15)

# add a legend
ax.legend(loc='best', frameon=False)

plt.savefig('graphics/exercise-2-1.png')
plt.show()