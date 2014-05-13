# suppose there is a 50% increase in delta
plt.figure(figsize=(8,6))
model.plot_phase_diagram(gridmax=100, N=1000, arrows=True, param='delta', 
                         shock=1.5, reset=True)
plt.show()