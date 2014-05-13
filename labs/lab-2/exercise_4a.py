# suppose there is a 100% increase in 1 / theta
plt.figure(figsize=(8,6))
model.plot_phase_diagram(gridmax=100, N=1000, arrows=True, param='theta', 
                         shock=0.5, reset=True)
plt.show()