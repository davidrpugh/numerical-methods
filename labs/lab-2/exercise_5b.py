model.plot_impulse_response(variables='all', method='linearization',
                            kind='efficiency_units', param='g', shock=0.5, T=50,
                            figsize=(8,12), log=False)
plt.show()