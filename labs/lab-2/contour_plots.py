from __future__ import division
import numpy as np
from scipy import optimize
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def F(K, L, A=1):
    """
    
    Classic constant elasticity of substitution (CES) production function.  
    Output is a function of the capital stock, K, and labor, L.  Note that 
    technology, A, is assumed to be labor-augmenting.
        
    """
    rho = (sigma - 1) / sigma
    
    # nest Cobb-Douglas as special case 
    if rho == 0:
        Y = K**alpha * (A * L)**(1 - alpha)
    else:
        Y =  (alpha * K**rho + (1 - alpha) * (A * L)**rho)**(1 / rho)
        
    return Y

def real_wage(K, L, A=1):
    """

    Assuming perfect competition, factors are paid there marginal products. Thus
    the real wage paid to labor is dF/dL.

    """
    rho = (sigma - 1) / sigma
    
    # nest Cobb-Douglas as special case
    if rho == 0:
        MPL = (1 - alpha) * (K / (A * L))**alpha
    else:
        MPL =  ((1 - alpha) * (A * L)**(rho - 1) * 
                (alpha * K**rho + (1 - alpha) * (A * L)**rho)**((1 / rho) - 1))
        
    return MPL

def return_capital(K, L, A=1):
    """

    Assuming perfect competition, factors are paid there marginal products. Thus
    the real return to capital is dF/dK less depreciation.

    """
    rho = (sigma - 1) / sigma
    
    # nest Cobb-Douglas as special case
    if rho == 0:
        MPK = alpha * (K / (A * L))**(alpha - 1)
    else:
        MPK =  (alpha * K**(rho - 1) * 
                (alpha * K**rho + (1 - alpha) * (A * L)**rho)**((1 / rho) - 1))
        
    return MPK - delta

def expansion_path(K, A):
    """

    Set of tangency points between isocost lines and isoquants. 

    """
    path = lambda L: (r / w) - (return_capital(K, L, A) / real_wage(K, L, A))
    
    res = optimize.root(path, x0=K, method='hybr')
    
    return res.x
    
def C(K, L):
    """

    Two sources of costs: labor costs, wL, and rental costs of capital, rK. 

    """
    return w * L + (r + delta) * K

def Pi(K, L, A):
    """

    Profits function from firm's maximization problem.

    """
    return F(K, L, A) - C(K, L)

# create a grid of (x,y) values which we will pass to function
capital = np.linspace(0, 50, 100)
labor   = np.linspace(0, 50, 100)
K, L    = np.meshgrid(capital, labor)

# Choose parameter values
alpha = 0.33
sigma = 1.0
delta = 0.04
A     = 1.0
w, r  = real_wage(2, 1, A), return_capital(2, 1, A)

# we will eventually plot...
output  = F(K, L, A)
costs   = C(K, L)
profits = Pi(K, L, A)

##### Contour plots for output and costs #####

# create a new Figure object 
fig = plt.figure(figsize=(12,6))

# create the contour plot for the second subplot
ax1 = fig.add_subplot(121)
production_surface = ax1.imshow(output, interpolation='gaussian', 
                                origin='lower', cmap=mpl.cm.hot, 
                                vmin=0, vmax=np.max(output), 
                                extent=(0, 50, 0, 50))

# demarcate the contours
CS = ax1.contour(K, L, output, np.linspace(0, 50, 10), colors=np.repeat('k', 10), 
                 linewidths=1, linestyles='solid')
ax1.clabel(CS, inline=1, fmt='%1.2f')

# axes, labels, title, colorbar etc.
ax1.set_xlim(0, 50)
ax1.set_ylim(0, 50)
ax1.set_xlabel(r'Capital, $K$', family='serif', fontsize=15)
ax1.set_ylabel(r'Labor, $L$', family='serif', fontsize=15)
ax1.set_title('Isoquants', family='serif', fontsize=15)

# create the contour plot for cost function
ax2 = fig.add_subplot(122)
cost_surface = ax2.imshow(costs, interpolation='gaussian', origin='lower', 
                          cmap=mpl.cm.cool, vmin=0, vmax=np.max(costs), 
                          extent=(0, 50, 0, 50))

# demarcate the contours
CS = ax2.contour(K, L, costs, np.linspace(0, 50, 10), colors=np.repeat('k', 10), 
                 linewidths=1, linestyles='solid')
ax2.clabel(CS, inline=1, fmt='%1.2f')

# axes, labels, title, colorbar etc.
ax2.set_xlim(0, 50)
ax2.set_ylim(0, 50)
ax2.set_xlabel(r'Capital, $K$', family='serif', fontsize=15)
ax2.set_ylabel(r'Labor, $L$', family='serif', fontsize=15)
ax2.set_title('Isocosts', family='serif', fontsize=15)

fig.text(0.5, 0.925, 'Constant returns CES production with $\sigma=%g$' % sigma,
         fontsize=20, ha='center', family='serif')

plt.savefig('graphics/contour-plots-for-output-and-costs.png')
plt.show()

##### Firm Size is indeterminate with CRTS #####

# create a new Figure object 
fig = plt.figure(figsize=(8,8))

# create the contour plot for the second subplot
ax1 = fig.add_subplot(111)

# demarcate the contours
CS_output = ax1.contour(K, L, output, np.linspace(0, 50, 10), cmap=mpl.cm.hot,
                        vmin=0, vmax=np.max(output), linewidths=1, 
                        linestyles='solid')
ax1.clabel(CS_output, inline=1, fmt='%1.2f')

CS_costs = ax1.contour(K, L, costs, np.linspace(0, 50, 10), cmap=mpl.cm.cool,
                       vmin=0, vmax=np.max(costs), linewidths=1, 
                       linestyles='solid')

# add the expansion path
ax1.plot(capital, np.array([expansion_path(k, A)[0] for k in capital]), 'r', 
         label=r'$\frac{K}{L}$')
         
# axes, labels, title, colorbar etc.
ax1.set_xlim(0, 50)
ax1.set_ylim(0, 50)
ax1.set_xlabel(r'Capital, $K$', family='serif', fontsize=15)
ax1.set_ylabel(r'Labor, $L$', family='serif', fontsize=15)
ax1.set_title(r'With CRTS, firm size is indeterminate!', family='serif', 
              fontsize=20)
ax1.legend(loc='upper right', bbox_to_anchor=(1.15, 1.0), frameon=False)

plt.savefig('graphics/firm-size-indeterminate-with-CRTS.png')
plt.show()

##### Profits surface and contour plot #####

# create a new Figure object 
fig = plt.figure(figsize=(12,6))

# create a 3D Axes object for the first subplot
ax = fig.add_subplot(121, projection='3d', elev=30, azim=-120)

# note the use of the new plot command!
profits_surface = ax.plot_surface(K, L, profits, rstride=1, cstride=1, 
                                  cmap=mpl.cm.winter_r, linewidth=0, 
                                  vmin=np.min(profits), vmax=0, antialiased=False)

# add the expansion path (i.e., the zero contour!) 
zero_contour = ax.plot(capital, 
                        np.array([expansion_path(k, A)[0] for k in capital]),
                        np.zeros(capital.size),
                        'r', label=r'$\frac{K}{L}$')
                        
# axes, labels, title, colorbar etc.
ax.set_xlim(0, 50)
ax.set_ylim(0, 50)
ax.set_xlabel(r'Capital, $K$')
ax.set_ylabel(r'Labor, $L$')
ax.set_zlabel(r'Profits, $\Pi$')

# create the contour plot for the second subplot
ax2 = fig.add_subplot(122)
im = ax2.imshow(profits, interpolation='gaussian', origin='lower', cmap=mpl.cm.winter_r, 
                vmin=np.min(profits), vmax=0, extent=(0, 50, 0, 50))

# demarcate the contours
CS = ax2.contour(K, L, profits, np.linspace(np.min(profits), 0, 10), colors=np.repeat('k', 10), 
                 linewidths=1, linestyles='solid')
ax2.clabel(CS, inline=1, fmt='%1.2f')

# add the expansion path (i.e., the zero contour!) 
zero_contour = ax2.plot(capital, 
                        np.array([expansion_path(k, A)[0] for k in capital]),
                        'r', label=r'$\frac{K}{L}$')

# axes, labels, title, colorbar etc.
ax2.set_xlim(0, 50)
ax2.set_ylim(0, 50)
ax2.set_xlabel(r'Capital, $K$')
ax2.set_ylabel(r'Labor, $L$')

# add a color bar to the figure
fig.colorbar(profits_surface, shrink=0.75, aspect=10)

# add a title to the figure
fig.text(0.5, 0.95, 'Firms earn zero profits in equilibrium', fontsize=20, ha='center')

# save the plot!
#plt.savefig('2012-12-22-Zero-profits-with-CRTS.png')

# display the plot!
plt.show()
