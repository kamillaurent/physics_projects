import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

# Read data from the .dat file
data = np.loadtxt('charged_jets_q.dat')

# Extract columns
eta_values = data[:, 1]
phi_values = data[:, 2]
pt_values = data[:, 3]

# Define range and normalization scale
x_range = (-0.4, 0.4)
y_range = (-0.4, 0.4)
pt_min = 0
pt_max = 50

n_bins = 33
hist, x_edges, y_edges = np.histogram2d(
    eta_values,
    phi_values, 
    bins=n_bins, range=[x_range, y_range], weights=pt_values/1129.)

#standardize the pixel
#zero_hist = hist - hist.mean()

# Create estandardized_hist = zero_hist / (hist.std() + 1e-5)
#xtent for the plot
extent = [-0.4, 0.4, -0.4, 0.4]

# Remove axis indices
plt.xticks([])
plt.yticks([])

# Plot the histogram
plt.figure(figsize=(10, 10))  # Set the figure size in inches
plt.imshow(hist.T, extent=extent, origin='lower', aspect='auto', cmap='viridis')
plt.colorbar(label='pT')
plt.xlabel('translated Eta')
plt.ylabel('translated Phi')
plt.title('Mean pt per Pixel in Quark Initiated Jets')
plt.tight_layout()

# Save the plot as a JPG image
plt.savefig('mean_q_charged.jpg', dpi=100)  # dpi controls the resolution

# Show the plot (optional)
plt.show()