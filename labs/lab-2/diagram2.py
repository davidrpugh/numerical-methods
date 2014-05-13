# compute the stable manifold using forward shooting
k0 =  model.steady_state.values['k_star'] / 10
ms_lower = model.solve_forward_shooting(k0, h=1.0, tol=1e-4, integrator='dopri5')

k0 = 3.5 * model.steady_state.values['k_star']
ms_upper = model.solve_forward_shooting(k0, h=1.0, tol=1e-4, integrator='dopri5')

plt.figure(figsize=(8,6))

# plot the phase diagram
model.plot_phase_diagram(gridmax=130, N=1000, arrows=True)

# plot the stable manifold
model.plot_trajectory(ms_upper, color='r')
model.plot_trajectory(ms_lower, color='r')

# demarcate the initial condition
plt.axvline(k0, linestyle='dashed', color='k')
plt.xticks([k0], ['$k_0$'])

# set of initial conditions for c
N = 20
c_lower   = 0.99 * ms_upper[0,2]
c_upper   = 1.01 * ms_upper[0,2]
init_vals = np.linspace(c_lower, c_upper, N)

# color scheme
color_map = mpl.cm.cool(np.linspace(0, 1, N))

for i, c0 in enumerate(init_vals):
    
    # simulate the model
    traj = model.integrate(t0=0, y0=[k0, c0], h=0.1, T=300, integrator='dopri5')

    # plot the trajectory 
    model.plot_trajectory(traj, color=color_map[i])

# compute the unstable manifold
eps = 1e-5
k_star, c_star = model.steady_state.values['k_star'], model.steady_state.values['c_star']

mu_lower = model.integrate(t0=0, y0=[k_star + eps, c_star], h=0.1, T=300, integrator='dopri5')
mu_upper = model.integrate(t0=0, y0=[k_star, c_star + eps], h=0.1, T=300, integrator='dopri5')

# plot the unstable manifold
model.plot_trajectory(mu_upper, color='r')
model.plot_trajectory(mu_lower, color='r')

# change the plot title
plt.title('Deviations between $M_S (M_U)$ are magnified (squashed)!', fontsize=20, family='serif')
plt.savefig('graphics/deviations-from-MS-MU.png')
plt.show()