from __future__ import division
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

# create a loaded coin
loaded_coin = stats.distributions.bernoulli(0.75)

# flip the loaded coin
for T in [10, 100, 1000, 10000]:
    
    # create a loaded coin
    loaded_coin_flips = loaded_coin.rvs((T,))

    # count the fraction of heads...
    fraction_heads = np.mean(loaded_coin_flips)
    
    print 'Fraction of heads with T=%i trials is %.3f' %(T, fraction_heads)