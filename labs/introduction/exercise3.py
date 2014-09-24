import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

#Set degrees of freedom
df = 1

# create a t-distribution object
students_t_rv = stats.distributions.t(df)

# set the seed
np.random.seed(12)

# generate data
T = 1000
sample_path = students_t_rv.rvs((T,))

# create new Figure and Axes objects
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)

# plot the sample path
ax.plot(sample_path, 'g-')

# label axes
xlab = ax.set_xlabel('Index', fontsize=15)
xlab.set_family('serif')

ylab = ax.set_ylabel('Observation', fontsize=15)
ylab.set_family('serif')

# set the title
title = ax.set_title("Random sample of data from Student's t with %g d.f." % df,
                     fontsize=15)
title.set_family('serif')

plt.savefig('graphics/exercise-2-2.png')
plt.show()