import matplotlib.pyplot as plt

plt.figure(figsize=(10,8))
ax = model.plot_phase_diagram(70, cmap='winter', arrows=True)[0]

kstar = model.steady_state.values['k_star']
cstar = model.steady_state.values['c_star']

# plot the unstable manifold
init = [kstar - 1e-5, cstar + 1e-5]
traj = model.integrate(t0=0, y0=init, T=1e3, h=0.5, integrator='lsoda')
ax.plot(traj[:,1], traj[:,2], 'r--', label='$M_U$',)

ax.hlines(traj[::10,2],traj[::10,1], traj[::10,1] + 0.1 * kstar, colors='r')
ax.vlines(traj[::10,1],traj[::10,2], traj[::10,2] - 0.035 * cstar, colors='r')

init = [kstar + 1e-5, cstar - 1e-5]
traj = model.integrate(t0=0, y0=init, T=1e3, h=0.5, integrator='lsoda')
ax.plot(traj[:,1], traj[:,2], 'r--')

ax.hlines(traj[::10,2],traj[::10,1], traj[::10,1] - 0.1 * kstar, colors='r')
ax.vlines(traj[::10,1],traj[::10,2], traj[::10,2] + 0.035 * cstar, colors='r')

# plot the stable manifold
optimal_traj = model.solve_forward_shooting(0.1 * kstar, h=1.0, tol=1e-5)
ax.plot(optimal_traj[:,1], optimal_traj[:,2], 'r', label='$M_S$')

ax.hlines(optimal_traj[::10,2],optimal_traj[::10,1] - 0.1 * kstar, 
          optimal_traj[::10,1], colors='r')
ax.vlines(optimal_traj[::10,1],optimal_traj[::10,2] - 0.035 * cstar, 
          optimal_traj[::10,2], colors='r')

optimal_traj = model.solve_forward_shooting(6.0 * kstar, h=1.0, tol=1e-5)
ax.plot(optimal_traj[:,1], optimal_traj[:,2], 'r')

ax.hlines(optimal_traj[::10,2],optimal_traj[::10,1], 
          optimal_traj[::10,1] + 0.1 * kstar, colors='r')
ax.vlines(optimal_traj[::10,1],optimal_traj[::10,2], 
          optimal_traj[::10,2] + 0.035 * cstar, colors='r')

ax.set_title('Optimal growth model has stable and unstable manifolds', 
             fontsize=20)
ax.legend(loc='best', frameon=False, prop={'family':'serif'})

plt.savefig('graphics/ramsey-phase-diagram-with-manifolds.png')
plt.show()