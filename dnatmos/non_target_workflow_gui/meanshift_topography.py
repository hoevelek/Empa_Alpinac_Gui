from pathlib import Path
import numpy as np
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn.datasets import make_blobs
from alpinac_sups.read_h5 import read_desc, read_ionisation_modes, read_mass_axis, H5_Reader, read_time_axis


# Load the hillshade topography image into a 2D array
#topography = np.load("topography.npy")
#topography = np.load(r"C:\Users\kaho\Desktop\data\data_Empa\210623.1218.std.9.2Dpic_max_greater_than_neigbour_diag.jpeg")




file = Path(r"C:\Users\kaho\TOFWerk_data\shared_data\ExtractTofSpectra_example\lib\DataFile_2021.09.26-16h56m06s.h5")
file = Path(r"G:\503_Themen\Klima\DNAtmos\test_data221115\DataFile_2021.09.26-16h56m06s.h5")

ei_seg = (1,5)
ci_seg = (6,9)

file = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230401.1747.air.3.h5")
ei_seg = (3,4)
ci_seg = (1,2)


#file = Path(r"C:\Users\kaho\data_Empa\210623.1218.std.9.h5")
#ei_seg = (1,3)
#ci_seg = (3,5)
file_b = str(file).encode()

ionisation_dict = {
    # Last value is not included
    "EI": ei_seg,
    "CI": ci_seg,
}




# %% Reader class

reader = H5_Reader(

    file,
    ionisation_dict

)

# %% This will read only the EI data
data_EI = reader.read_full_ionization_data('EI')
# %% This will return the whole data over all segements
#data_all = reader.read_full_tof_data()
data = data_EI
#print data dimensions
print(data.shape)

# select subarray




topography = data[600:750, 600:750]

print(topography.shape)
#plot topography 3D
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x = np.arange(0, topography.shape[0], 1)
y = np.arange(0, topography.shape[1], 1)
X, Y = np.meshgrid(x, y)
ax.plot_surface(X, Y, topography)
plt.show()

# desity plot of topography
import seaborn as sns
sns.heatmap(topography)

#plot topography 2D
import matplotlib.pyplot as plt
plt.imshow(topography)




# # Define the mean shift window size
# window_size = 50

# # Apply mean shift to the image
# ms = MeanShift(window_size)
# ms.fit(topography)
# labels = ms.labels_
# cluster_centers = ms.cluster_centers_

# # Identify the cluster centers, which correspond to the hill tops
# hill_tops = []
# for center in cluster_centers:
#     row, col = np.where(topography == center[0])
#     hill_tops.append((row[0], col[0]))

# # Optionally, apply a threshold to filter out smaller hills
# threshold = 500
# filtered_hill_tops = [ht for ht in hill_tops if topography[ht] > threshold]

# # Plot the hill tops on the original image
# import matplotlib.pyplot as plt
# plt.imshow(topography, cmap='gray')
# plt.scatter([ht[1] for ht in filtered_hill_tops], [ht[0] for ht in filtered_hill_tops], s=50, c='r', marker='x')
# plt.show()


import numpy as np
from sklearn.cluster import MeanShift

# Load the hillshade topography image into a 2D array
#topography = np.load("topography.npy")

# Define the mean shift bandwidth
bandwidth = 50

# Apply mean shift to the image
topography_reshaped = topography.reshape(-1, 1)
bandwidth = estimate_bandwidth(topography_reshaped, quantile=0.2, n_samples=500)

ms = MeanShift(bandwidth=bandwidth)
ms.fit(topography_reshaped)
labels = ms.labels_
print(np.unique(labels))
cluster_centers = ms.cluster_centers_
n_clusters_ = len(np.unique(labels))



#get standard deviation from topography
std = np.std(topography)
mean_topography = np.mean(topography)
max_topography = np.max(topography)
threshold = std

import numpy as np

# Example 2D data of summed counts
summed_counts = np.array([[0, 0, 1, 0],
                          [2, 0, 0, 0],
                          [0, 0, 3, 1],
                          [0, 1, 0, 0]])

# Threshold value
threshold = 0.1

# Create a boolean mask
mask = topography > threshold

# Use the mask to select the corresponding elements
x_coords, y_coords = np.where(mask)
counts = topography[x_coords, y_coords]*1000000

# Create the output array
single_events_array = np.column_stack((x_coords, y_coords, counts))








# EXAMPLAE FROM SKLEARN

import matplotlib.pyplot as plt

#centers = [[1, 1], [-1, -1], [1, -1]]
#X, _ = make_blobs(n_samples=10000, centers=centers, cluster_std=0.6)
X = single_events_array

# The following bandwidth can be automatically detected using
bandwidth = estimate_bandwidth(X, quantile=0.2, n_samples=500)

ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
ms.fit(X)
labels = ms.labels_
cluster_centers = ms.cluster_centers_


labels_unique = np.unique(labels)
n_clusters_ = len(labels_unique)


plt.figure(1)
plt.clf()

# rainbow colors array with length of number of clusters
colors = plt.cm.rainbow(np.linspace(0, 1, n_clusters_))
markers = ['x']*n_clusters_

for k, col in zip(range(n_clusters_), colors):
    my_members = labels == k
    cluster_center = cluster_centers[k]
    plt.plot(X[my_members, 0], X[my_members, 1], markers[k], color=col)
    plt.plot(
        cluster_center[0],
        cluster_center[1],
        markers[k],
        markerfacecolor=col,
        markeredgecolor="k",
        markersize=14,
    )
plt.title("Estimated number of clusters: %d" % n_clusters_)
plt.show()








# plot topography including cluster centers
import matplotlib.pyplot as plt
plt.imshow(topography, cmap='gray')
plt.scatter(cluster_centers, [0 for i in range(len(cluster_centers))], s=50, c='r', marker='x') 
plt.show()


# import numpy as np
# from sklearn.cluster import MeanShift
# import matplotlib.pyplot as plt
# from PIL import Image

# # Load the grayscale image
# img = np.array(Image.open('path/to/image.png').convert('L'))

# # Reshape the image into a 2D array
# data = img.reshape((-1, 1))

# # Determine cluster centers using meanshift
# ms = MeanShift()
# ms.fit(data)
# labels = ms.labels_
# centers = ms.cluster_centers_

# Convert labels to 2D array for plotting
labels_2d = labels.reshape(topography.shape)
labels_flat = labels_2d.flatten()

#plot labels 2d
plt.imshow(labels_2d, cmap='gray')

# Create scatter plot of image data with colored points
fig, ax = plt.subplots()
#ax.scatter(range(topography.shape[1]), range(topography.shape[0]), c=labels_flat, cmap='viridis', alpha=0.5)
ax.imshow(topography, cmap='gray')

# Overlay cluster centers on scatter plot
for center in cluster_centers:
    ax.plot(center[0], center[1], 'k+', markersize=10, c= 'r')

# Show the plot
plt.show()












# Identify the cluster centers, which correspond to the hill tops



hill_tops = []
for center in cluster_centers:
    row, col = np.where(topography == center[0])
    hill_tops.append((row[0], col[0]))

# Optionally, apply a threshold to filter out smaller hills
threshold = 500
filtered_hill_tops = [ht for ht in hill_tops if topography[ht] > threshold]

# Plot the hill tops on the original image
import matplotlib.pyplot as plt
plt.imshow(topography, cmap='gray')
plt.scatter([ht[1] for ht in filtered_hill_tops], [ht[0] for ht in filtered_hill_tops], s=50, c='r', marker='x')
plt.show()
