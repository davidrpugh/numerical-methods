# suppose there is a 10% decrease in sigma
plt.figure(figsize=(8,6))
model.plot_phase_diagram(gridmax=50, N=1000, arrows=True, param='sigma', 
                         shock=0.9, reset=True)
plt.show()