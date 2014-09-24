# create new Figure object
fig = plt.figure(figsize(10, 12))

# create the first subplot
ax1 = fig.add_subplot(321)

# add the histogram of the first row of sample averages.
n, bins, patches = ax1.hist(sample_averages[0,:], bins=np.arange(0.5, 7, 1), 
                            align='mid', normed=True)

# don't forget to label axes!
ylabel = ax1.set_ylabel('Density')
ylabel.set_family('serif')

# add a title!
title = ax1.set_title('T=1')
title.set_family('serif')

# create the second subplot
ax2 = fig.add_subplot(322)

n, bins, patches = ax2.hist(sample_averages[49,:], bins=25, align='mid', normed=True)

x = np.linspace(0, 7, 1000)
normal = stats.distributions.norm(loc=mu, scale=sigma / np.sqrt(49))
ax2.plot(x, normal.pdf(x), 'r--', label=r'$N\left(\mu,\ \frac{\sigma}{\sqrt{T}}\right)$')

# don't forget to label axes!
ylabel = ax2.set_ylabel('Density')
ylabel.set_family('serif')

# set the title
ax2.set_title('T=50')

# create the third subplot
ax3 = fig.add_subplot(323)

n, bins, patches = ax3.hist(sample_averages[1,:], bins=np.arange(0.5, 7, 1), align='mid', normed=True)

x = np.linspace(0, 7, 1000)
normal = stats.distributions.norm(loc=mu, scale=sigma / np.sqrt(2))
ax3.plot(x, normal.pdf(x), 'r--', label=r'$N\left(\mu,\ \frac{\sigma}{\sqrt{T}}\right)$')

# don't forget to label axes!
ylabel = ax3.set_ylabel('Density')
ylabel.set_family('serif')

# set the title
ax3.set_title('T=2')

# create the fourth subplot
ax4 = fig.add_subplot(324, sharex=ax1)

n, bins, patches = ax4.hist(sample_averages[499,:], bins=25, align='mid', normed=True)

x = np.linspace(0, 7, 1000)
normal = stats.distributions.norm(loc=mu, scale=sigma / np.sqrt(500))
ax4.plot(x, normal.pdf(x), 'r--', label=r'$N\left(\mu,\ \frac{\sigma}{\sqrt{T}}\right)$')

# don't forget to label axes!
ylabel = ax4.set_ylabel('Density')
ylabel.set_family('serif')

# set the title
ax4.set_title('T=500')

# create the fifth subplot
ax5 = fig.add_subplot(325, sharex=ax1)

n, bins, patches = ax5.hist(sample_averages[2,:], bins=np.arange(0.5, 7, 1), align='mid', normed=True)

x = np.linspace(0, 7, 1000)
normal = stats.distributions.norm(loc=mu, scale=sigma / np.sqrt(3))
ax5.plot(x, normal.pdf(x), 'r--', label=r'$N\left(\mu,\ \frac{\sigma}{\sqrt{T}}\right)$')

# don't forget to label axes!
ax5.set_xlabel(r'$\hat{\mu}$')
ylabel = ax5.set_ylabel('Density')
ylabel.set_family('serif')

# set the title
ax5.set_title('T=3')

# create the sixth subplot
ax6 = fig.add_subplot(326, sharex=ax1)

n, bins, patches = ax6.hist(sample_averages[999,:], bins=25, align='mid', normed=True)

x = np.linspace(0, 7, 1000)
normal = stats.distributions.norm(loc=mu, scale=sigma / np.sqrt(1000))
ax6.plot(x, normal.pdf(x), 'r--', label=r'$N\left(\mu,\ \frac{\sigma}{\sqrt{T}}\right)$')

# don't forget to label axes!
ax6.set_xlabel(r'$\hat{\mu}$')
ylabel = ax6.set_ylabel('Density')
ylabel.set_family('serif')

# set the title
ax6.set_title('T=1000')

# add a title for the entire set of plots!
plt.figtext(0.5, 1.025, 'CLT in Action!', ha='center', fontsize=20, 
            family='serif')

# set the layout
plt.tight_layout()

plt.savefig('graphics/The-CLT-in-action.png')
plt.show()