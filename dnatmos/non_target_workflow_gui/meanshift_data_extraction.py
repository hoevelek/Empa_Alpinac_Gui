#%%
#import sys
#sys.path.append('c:\\Users\\kaho\\polybox\\ALPINAC\\Alpinac_packing\\Alpinac_packing\\alpinac_sups')
import time
t_start = time.time()
from pathlib import Path

from alpinac_sups.read_h5 import read_desc, read_ionisation_modes, read_mass_axis, H5_Reader, read_time_axis
from matplotlib.lines import segment_hits
import numpy as np
import matplotlib.pyplot as plt


#
#file = Path(r"C:\Users\coli\Documents\TW_EI_CI_data\DataFile_2021.09.26-16h56m06s.h5")


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


#%% specify an ionization idctionary (our files do not have one)

ionisation_dict = {

    # Last value is not included

    "EI": (1,5),

    "CI": (6,10),

}

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

data_all = reader.read_full_tof_data()


data = data_EI
data[0,1]

type(data)



#%%
data = data_EI

print("coordinates of maximal values in data")
max_coordinates = np.unravel_index(np.argmax(data), data.shape)
print(max_coordinates)

#%%

# Subset of data +-50 around the maximum
data_sub = data[max_coordinates[0]-20:max_coordinates[0]+20,max_coordinates[1]-50:max_coordinates[1]+50]
# print data_sub dimensions
data_sub2 = data[max_coordinates[0]-2:max_coordinates[0]+2,max_coordinates[1]-2:max_coordinates[1]+2]

print(data_sub2.shape)
plt.imshow(data_sub)
data = data_sub

#%%
# Prepare data for feature input
# load black and white image as grayscale
#image = cv2.imread('black_and_white_image.png', cv2.IMREAD_GRAYSCALE)

image = data

# get image dimensions
height, width = image.shape

# create x and y value arrays
x_values = np.arange(width)
y_values = np.arange(height)

# create meshgrid of x and y values
xx, yy = np.meshgrid(x_values, y_values)

# flatten image matrix to one-dimensional array
intensities = image.flatten()

# create feature matrix
feature_matrix = np.column_stack((xx.flatten(), yy.flatten(), intensities))



# print feature matrix
print(feature_matrix)



#%%
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn.datasets import make_blobs
# Meanshift algorithm on black and white image to find two clusters
# https://scikit-learn.org/stable/auto_examples/cluster/plot_mean_shift.html#sphx-glr-auto-examples-cluster-plot-mean-shift-py

# Generate sample data
#centers = [[1, 1], [-1, -1], [1, -1]]
#X, _ = make_blobs(n_samples=10000, centers=centers, cluster_std=0.6)
# prepare data for clustering
#data = X[:, ::-1]  # flip axes for better plotting
# convert array of shape [x,y] to [x1,y1,z1]

# Converting image into array of dimension [nb of pixels in originImage, 3]


#use only data points with intensity > 1% of the maximum

data_for_clustering = feature_matrix[feature_matrix[:,2] > 0.1*np.max(feature_matrix[:,2])]

#data_for_clustering = feature_matrix
print(data_for_clustering.shape)
print(data.shape)

# Compute clustering with MeanShift
# The following bandwidth can be automatically detected using
bandwidth = estimate_bandwidth(data_for_clustering, quantile=0.3, n_samples=1000)
print(bandwidth)
ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
ms.fit(data_for_clustering)
labels = ms.labels_
cluster_centers = ms.cluster_centers_
print(cluster_centers)

labels_unique = np.unique(labels)
n_clusters_ = len(labels_unique)

print("number of estimated clusters : %d" % n_clusters_)
print("cluster centers: ", cluster_centers)


# get cluster centers and their intensity values
#cluster_centers = ms.cluster_centers_
intensities = cluster_centers[:, 2]

# filter cluster centers by intensity value to find maxima
maxima = cluster_centers[intensities.argmax()]

# print coordinates of maxima
print("Maxima coordinates: ({}, {})".format(int(maxima[0]), int(maxima[1])))

#%%
# plot subset and cluster centers
import matplotlib.pyplot as plt
plt.imshow(data)
plt.scatter(cluster_centers[:,0], cluster_centers[:,1], c='r')
plt.show()

#plot projection onto x axis
plt.plot(data[:,int(maxima[1])])
#plot line at maxima x position
plt.axvline(x=int(maxima[1]), color='r')
plt.show()

#plot projection onto y axis
plt.plot(data[int(maxima[0]),:])
#plot line at maxima y position
plt.axvline(x=int(maxima[0]), color='r')




#%%
# COSINE SIMILARITY
import sklearn.metrics as metrics
metrics.pairwise.cosine_similarity(data, data)

# ALTERNATIVE PACKAGES

# plot subset as 2D color plot
import matplotlib.pyplot as plt
plt.imshow(data)
plt.show()

#%%
import sys
sys.path.append(r"C:\Users\kaho\polybox\DNAtmos_code\second_party_code\MeanShift_py")


import mean_shift as ms

#%%

#data = get_data_from_somewhere()
mean_shifter = ms.MeanShift()
mean_shift_result = mean_shifter.cluster(data, kernel_bandwidth = 1000)
original_points =  mean_shift_result.original_points
shifted_points = mean_shift_result.shifted_points
cluster_assignments = mean_shift_result.cluster_ids

print(original_points)
print(cluster_assignments)

#%%

# If you want to use multivariate gaussian kernel
# By default it uses unviariate gaussian kernel
# Make sure the dimensions of 'data' and the kernel match
mean_shifter = ms.MeanShift(kernel='multivariate_gaussian')
mean_shift_result = mean_shifter.cluster(data, kernel_bandwidth = [10,20,30])

print(mean_shift_result.shifted_points)

# %%


import mean_shift as ms
import matplotlib.pyplot as plt
import numpy as np

#data = np.genfromtxt('data.csv', delimiter=',')

mean_shifter = ms.MeanShift()
mean_shift_result = mean_shifter.cluster(data, kernel_bandwidth = 1)

original_points =  mean_shift_result.original_points
shifted_points = mean_shift_result.shifted_points
cluster_assignments = mean_shift_result.cluster_ids

#%%
x = original_points[:,0]
y = original_points[:,1]
Cluster = cluster_assignments
centers = shifted_points

fig = plt.figure()
ax = fig.add_subplot(111)
scatter = ax.scatter(x,y,c=Cluster,s=50)
for i,j in centers:
    ax.scatter(i,j,s=50,c='red',marker='+')
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.colorbar(scatter)

fig.savefig("mean_shift_result")

# %%
