import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

# Read data from the .dat file
data_path = 'example1.dat'
data = open(data_path, 'r')

# Initialize variables
eta_values = []
phi_values = []
pt_values = []
group_count = 0

# Process each line in the file
for line in data:
    line = line.strip()  # Remove leading/trailing whitespace
    if line:  # Non-empty line
        # Split the line into values
        values = line.split()
        eta_values.append(float(values[1]))
        phi_values.append(float(values[2]))
        pt_values.append(float(values[3]))       
    else:  # Blank line, indicates the end of a group
        if group_count > 0:
            # Center eta and phi in the pT centroid
            #eta_centroid = np.sum(np.array(eta_values) * np.array(pt_values)) / np.sum(np.array(pt_values))
            #phi_centroid = np.sum(np.array(phi_values) * np.array(pt_values)) / np.sum(np.array(pt_values))

            #eta_values -= eta_centroid
            #phi_values -= phi_centroid

            # define ranges and colormap
            x_range = (-0.4, 0.4)
            y_range = (-0.4, 0.4)
            pt_min = 0
            pt_max = 100

            n_bins = 33
            hist, x_edges, y_edges = np.histogram2d(eta_values, phi_values, bins=n_bins, range=[x_range, y_range], weights=pt_values)

            # Create standardized histogram
            #zero_hist = hist - hist.mean()
            #standardized_hist = zero_hist / (hist.std() + 1e-5)

            # Create extent for the plot
            extent = [-0.4, 0.4, -0.4, 0.4]

            # Set the colormap range for the color normalization
            norm = Normalize(vmin=pt_min, vmax=pt_max)


            # Plot the histogram
            plt.figure(figsize=(33, 33))
            plt.imshow(
                hist.T, extent=extent, 
                origin='lower', 
                aspect='auto', cmap='viridis',
                norm = norm
            )
            #plt.colorbar(label='pT')
            #plt.xlabel('Eta')
            #plt.ylabel('Phi')
            #plt.title('2D Histogram of pT in Phi-Eta Space')
            #plt.tight_layout()

            # Remove axis indices and labels
            plt.xticks([])
            plt.yticks([])

            # Save the plot as a JPG image in the specified folder
            image_name = 'dataset/' + f'quark_jet_{group_count}.jpg'
            plt.savefig(image_name, dpi=100, bbox_inches='tight', pad_inches=0)

            # Close the figure
            plt.close()

            # Show the plot (optional)
            # plt.show()

            # Reset variables for the next group
            eta_values = []
            phi_values = []
            pt_values = []
        
        group_count += 1

# Close the file after processing
data.close()