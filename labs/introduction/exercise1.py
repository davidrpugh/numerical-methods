import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

# specify degrees of freedom
df = 3

# create a random variable object representing the normal distribution
chi_square_rv = stats.distributions.chi2(df)

# check type of chi_square_rv
type(chi_square_rv)

# create a grid of values at which to evaluate pdf and cdf
grid = np.linspace(0, 15, 1000)

# create new Figure object
fig = plt.figure(figsize=(8, 10))

# plot the pdf
ax = fig.add_subplot(211)

ax.plot(grid, chi_square_rv.pdf(grid), 'r-')
ax.axis('tight')
ax.set_ylim(0, 0.5)

# label the axes...
xlab = ax.set_xlabel('X', fontsize=15)
xlab.set_family('serif')
ylab = ax.set_ylabel('f(x)', fontsize=15, rotation='horizontal')
ylab.set_family('serif')

# ...and add a title
ax.set_title(r'Probability density function of a $\chi^2$ r.v. with %g d.f.' %df,
            family='serif', fontsize=15)

# plot the cdf
ax1 = fig.add_subplot(212)

ax1.plot(grid, chi_square_rv.cdf(grid), 'b-')
ax1.axis('tight')

# label the axes...
xlab = ax1.set_xlabel('X', fontsize=15)
xlab.set_family('serif')
ylab = ax1.set_ylabel('F(x)', fontsize=15, rotation='horizontal')
ylab.set_family('serif')

# ...and add a title
ax1.set_title(r'Distribution function of a $\chi^2$ r.v. with %g d.f.' %df,
              family='serif', fontsize=15)

# adjust the layout
fig.tight_layout()

plt.savefig('graphics/exercise-1-1a.png')
plt.show()

# create new Figure and Axes objects
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)

# plot the survival function
ax.plot(grid, chi_square_rv.sf(grid), 'g-')
ax.set_xscale('log')
ax.set_yscale('log')

# label the axes...
xlab = ax.set_xlabel('X (log scale)', fontsize=15)
xlab.set_family('serif')
ylab = ax.set_ylabel('Survival function, 1 - F(x) (log scale)', fontsize=15)
ylab.set_family('serif')

# ...and add a title
ax.set_title(r'Survival function of a $\chi^2$ r.v. with %g d.f.' %df, 
             family='serif', fontsize=15)

plt.savefig('graphics/exercise-1-1b.png')
plt.show()