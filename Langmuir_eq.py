import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define the Langmuir equation
def langmuir_equation(x, B_max, K_d):
    return (B_max * x) / (K_d + x)

start = 0.00
end = 300
increment = 15

array = [start + i * increment for i in range(int((end - start) / increment) + 1)]
print(array)

# Your experimental data (replace with your actual data)
concentrations = array  # Concentrations of target added at each titration
fluorescence_intensities = np.array([68.63, 110.18,	293.49,	485.19,	701.14,	914.75,	1145.79,	1381.48,	1616.81,	1866.87,	2074.75,	2307.47,	2496.07,	2670.10,	2785.10,	2898.67,	2997.86,	2994.60,	2993.51,	2985.81,	2990.24]) 
#Short MB ([68.63, 110.18,	293.49,	485.19,	701.14,	914.75,	1145.79,	1381.48,	1616.81,	1866.87,	2074.75,	2307.47,	2496.07,	2670.10,	2785.10,	2898.67,	2997.86,	2994.60,	2993.51,	2985.81,	2990.24])  # Corresponding fluorescence intensities

#([89.405,	94.776,	114.63,	175.274,	208.58,	266.532,	335.747,	392.301,	461.782,	484.551,	528.12,	532.692,	552.895,	587.444,	582.545,	578.224,	578.352,	562.278,	538.686,	506.203,	497.155	])

# Initial estimates for B_max and K_d
B_max_initial = max(fluorescence_intensities)
K_d_initial = 1  # Replace with your best estimate for K_d

# Fit the data to the Langmuir equation
popt, _ = curve_fit(langmuir_equation, concentrations, fluorescence_intensities, p0=[B_max_initial, K_d_initial])
B_max, K_d = popt

# Plot the experimental data
plt.scatter(concentrations, fluorescence_intensities, color='red', label='Experimental data')

# Plot the fitted Langmuir equation
x_values = np.linspace(min(concentrations), max(concentrations), 100)
y_values = langmuir_equation(x_values, B_max, K_d)
plt.plot(x_values, y_values, label=f'Langmuir Fit (B_max={B_max:.2f}, K_d={K_d:.2f})')

# Customize the plot
plt.xlabel('Concentration of target')
plt.ylabel('Fluorescence intensity')
plt.legend()
plt.grid(False)
plt.title('Fluorescence Spectroscopy Titration')

plt.savefig('Langmuir_LongTarget.svg', dpi=600, bbox_inches='tight')

# Show the plot
plt.show()

