# Script executed in virtual environment, change to it via:  PS c:/Users/kaho/polybox/DNAtmos_code/venv_dnatmos/Scripts/activate.ps1


# STEP 1: do extraction of data from h5 file and save as txt file
# see code below
# STEP 2: do mode identification 
# py -m alpinac.mode_identification C:\Users\kaho\polybox\ALPINAC\Alpinac_packing\halosearch\data\recetox_exposome\RECETOX_GC-EI-MS_20201028\1_2_3_7_8_Pentachlorodibenzo_p_dioxin__C12H3Cl5O2.txt --fitting-method discretepy -m alpinac.mode_identification C:\Users\kaho\polybox\ALPINAC\Alpinac_packing\halosearch\data\recetox_exposome\RECETOX_GC-EI-MS_20201028\1_2_3_7_8_Pentachlorodibenzo_p_dioxin__C12H3Cl5O2.txt --fitting-method discrete


#%% install packages from repository


from platform import python_version
from pathlib import Path
import sys
import time
import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.patches as patches
import matplotlib.lines as lines
from matplotlib.lines import segment_hits


#%%
from platform import python_version
print(python_version())
import sys
print(sys.path)




#%% install local packages
# e.g. via pip install -e C:\Users\kaho\polybox\ALPINAC\Alpinac_packing\Alpinac_packing\alpinac_sups
# or via pip install -e C:\Users\kaho\polybox\ALPINAC\Alpinac_packing\Alpinac_packing\halosearch
# from alpinac.install_packages import install_packages # TODO install all requierements if not alraedy installed by importing alpinac
# install_packages()
import alpinac


#%%
# sys.path.append('c:\\Users\\kaho\\polybox\\ALPINAC\\Alpinac_packing\\Alpinac_packing\\alpinac_sups') # only if not already in path
from alpinac_sups.read_h5 import H5_Reader, read_time_axis, read_desc, read_ionisation_modes, read_mass_axis



# import sys
import time
t_start = time.time()

#file = Path(r"C:\Users\coli\Documents\TW_EI_CI_data\DataFile_2021.09.26-16h56m06s.h5")


file = Path(r"C:\Users\kaho\Desktop\data\TOFWerk_data\shared_data\ExtractTofSpectra_example\lib\DataFile_2021.09.26-16h56m06s.h5")
#file = Path(r"G:\503_Themen\Klima\DNAtmos\tmp_data\test_data221115\DataFile_2021.09.26-16h56m06s.h5")
ei_seg = (1,5)
ci_seg = (6,9)


#file = Path(r"C:\Users\kaho\data_Empa\210623.1218.std.9.h5")
#ei_seg = (1,3)
#ci_seg = (3,5)
file_b = str(file).encode()




#ionisation_dict = {
#    # Last value is not included
#    "EI": ei_seg,
#    "CI": ci_seg,
#}

file = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230401.1747.air.3.h5")
ei_seg = (3,4)
ci_seg = (1,2)


#%% specify an ionization idctionary (our files do not have one)

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

data = reader.read_full_ionization_data('CI')

# %% This will return the whole data over all segements

data_all = reader.read_full_tof_data()

print(data.size)

print(type(data))

# %% Get mass calibration parameters
#Mass calibration parameters
masscalibration_p1 = 2829.4940436162415
masscalibration_p2 = -4042.695405928985
masscalibration_p3 = 0.49993765722959954

# TODO add mass calibration parameter to reader class
# print attributes of reader
print(reader.__dict__)
#reader.__getattribute__('MassCalibration_p1')








    ### HERE STARTS EXPERIMENTAL CODE TO FIND MAXIMA FROM IMAGE 


    # STRATEGY:
    # find gridd on which potentially masses and RT could have max in a summed 1D representation
    # to check if a retention time peak is superposition try to check if fits for all masses for this RT, define criteria
    # Potentially use BillerBiemannAlgo for RT, but only summed 1D!!
    # a transformation in a symmetric distribution could help
    
    # go to each crossing point of gridd and a) find region by going to neighbours and add them to gridd point as long as above threshold.
    # e.g. go away from center. Define diagonal and straight steps with two increments in horizontal and vertical direction
    # step(-1,-1): case diagonal go to three next pixels step(0,-1), step(-1, 0), step(-1,-1)
    # step(-1,0): case orthogonal go to lower next pixels step(-1,0)

    
    # To determining noise more properly, Lionels filter could also use diagonal neighbours (a new operation, for finding grid, this is better). This would
    # Do special treatment if coeluting. This will be the case if the next griddpoint is within the pixels that were added.
    # what this special treatment should be should be discussed, maybe shape analysis of border? Maybe fit? Maybe Lionels filter with small width?

    #Which parameter would we have?
    # 2x 2 from peakfinder
    # evtl from Billman alg.
    # LOD full intensity
    # LOD binary



    # idea: convert to picture, do image manipulation at GPU
    # e.g. use OpenCl, pyopencl





    #from matplotlib import cm
    # load tof data: 
    # y dimension: tof
    # x dimension: retention time?


    # FILTER TO DETECT VERTICL PLUS HORIZONTAL STEPS AND RETURN BINARY=1 IF DETECTED, 0 ELSE


#%% import packages 

import numpy as np
from PIL import Image, ImageOps, ImageFilter
from scipy.ndimage import maximum_filter, gaussian_filter


# FILTER TO DETECT VERTICAL PLUS HORIZONTAL STEPS AND RETURN BINARY=1 IF DETECTED, 0 ELSE

#%% import packages 


N_NEIGHBORS = 4
def find_max_2D(a):
    # Will store the wether of all the neighbors
    # of a peak will be larger than them
    greater_than_neighbors = np.ones_like(a, dtype=bool)
    for i in range(1, N_NEIGHBORS + 1):
        # x axis
        greaters_0 = np.roll(a, i, axis=0) < a
        greaters_1 = np.roll(a, -i, axis=0) < a
        # y axis
        greaters_2 = np.roll(a, i, axis=1) < a
        greaters_3 = np.roll(a, -i, axis=1) < a

        greater_than_neighbors = np.logical_and.reduce(
            [
                greater_than_neighbors,
                greaters_0,
                greaters_1,
                greaters_2,
                greaters_3,
            ],
            axis=0,
            dtype=bool,
        )
    return greater_than_neighbors

#%% do filtering

N_NEIGHBORS_DIAG = 4
def find_max_2D_diag(a):
# Will store the wether of all the neighbors
# of a peak will be larger than them
    greater_than_neighbors_diag = np.ones_like(a, dtype=bool)
    #for i in range(1, N_NEIGHBORS_DIAG + 1):
    for i in range(1, N_NEIGHBORS_DIAG + 1):
        # x axis, y axis: right, up
        greaters_0 = np.roll(np.roll(a, i, axis=0), i, axis=1) <= a
        # x,y: right, down
        greaters_1 = np.roll(np.roll(a, i, axis=0), -i, axis=1) <= a
        # x,y axis: left, up
        greaters_2 =  np.roll(np.roll(a, -i, axis=0), i, axis=1) <= a
        # x,y axis: left, down
        greaters_3 = np.roll(np.roll(a, -i, axis=0), -i, axis=1) <= a

        greater_than_neighbors_diag = np.logical_and.reduce(
            [
                greater_than_neighbors_diag,
                greaters_0,
                greaters_1,
                greaters_2,
                greaters_3,
            ],
            axis=0,
            dtype=bool,
        )
    return greater_than_neighbors_diag



#%% extrac peaks

# EXTRACT PEAKS

#convert to 8bit image

# Image has x,y,

im_raw = Image.fromarray(data).convert("L")

 # SMOOTH IN 2D
im_autocontrast_filter = im_raw.filter(ImageFilter.SMOOTH)#.filter(ImageFilter.FIND_EDGES)

im_autocontrast_filter.save(Path(r"C:\Users\kaho\Desktop\data\data_Empa\210623.1218.std.9.2Dpic_autocontrast_smoothed.jpeg"))


#%%

# CONVERT BACK TO NP ARRAY (REMOVE OTHER ARRAYS FROM MEM IF READY!)
filtered_im_as_np = np.asarray(im_autocontrast_filter)
print("shape image")
print(filtered_im_as_np.shape)
hori_dim, vert_dim = filtered_im_as_np.shape



# %%
# GENERATE RT AND MASS SUM SPECTRA OF INT8 IMAGE DATA - SMOOTHED FOR FULL INTENSITY

# find stripes, all is binary
vertical_lines = np.sum(filtered_im_as_np, axis = 0)
horizontal_lines = np.sum(filtered_im_as_np, axis = 1)
max_vert = np.max(vertical_lines)
max_hori = np.max(horizontal_lines)


from scipy.signal import find_peaks
vert_max = find_peaks(vertical_lines, prominence = (max_vert*0.0001, None), distance =7)
hori_max = find_peaks(horizontal_lines, prominence = (max_hori*0.0001, None), distance =7)


# %%
def plot_lines(vertical_lines, vert_max):
    import matplotlib.pyplot as plt
    plt.plot(vertical_lines)
    plt.plot(vert_max[0], vertical_lines[vert_max[0]], "x")
    plt.show()

fig, ax = plt.subplots()
plot_lines(vertical_lines, vert_max[0])
fig, ax = plt.subplots()
plot_lines(horizontal_lines, hori_max[0])

plot_lines(vertical_lines, vert_max)
plot_lines(horizontal_lines, hori_max)

#STOP

# %%
# GENERATE RT AND MASS SUM SPECTRA OF INT8 IMAGE DATA - SMOOTHED FOR NUMBER DENSITY SPECTRA

# binary filter enhanced picture
bin_enhanced = find_max_2D(filtered_im_as_np)
print(bin_enhanced.shape)

# find stripes, all is binary
vertical_lines = np.sum(bin_enhanced, axis = 0)
horizontal_lines = np.sum(bin_enhanced, axis = 1)
max_vert = np.max(vertical_lines)
max_hori = np.max(horizontal_lines)


from scipy.signal import find_peaks
vert_max = find_peaks(vertical_lines, prominence = (max_vert*0.015, None), distance =7)
hori_max = find_peaks(horizontal_lines, prominence = (max_hori*0.05, None), distance =7)

plt.figure()
plot_lines(vertical_lines, vert_max)
plt.title("vertical lines")
plt.figure()
plot_lines(horizontal_lines, hori_max)
plt.title("horizontal lines")


# %%
vert_max_pos = vert_max[0]
print(type(vert_max_pos))
vert_max_mass_regime = np.power(vert_max_pos, 2)+0
counts, bins = np.histogram(vert_max_mass_regime, bins = 100)

vert_max_prom = vert_max[1]['prominences']
print(vert_max_prom)
plt.figure()
plt.plot(vert_max_mass_regime,vert_max_prom)
plt.title("vertical lines, mass regime")


# %%
# AUTO-MASS CALIBRATION
# calculate_distances_of_entries_of_1D_array
print(vert_max)
# print prominance of peaks
print(vert_max[1]["prominences"])
# print positions of peaks
print(vert_max[0])



# %%
# calculate distances to n nearest neighbors

def calculate_distances_to_n_nearest_neighbors_of_1D_array(position_bin, prominences, n_neighbors):
    distances = []
    bins = []
    weights = []

    #array = array_2D[:,0]
    #bin_no = array_2D[:,1]
    for i in range(len(position_bin)):
        print(f"i={i}")
        for n in range(1, n_neighbors+1):
            print(n)
            if i+n < len(position_bin):
                distances.append(-position_bin[i]+position_bin[i+n])
                bins.append(0.5*(position_bin[i]+position_bin[i+n]))
                weights.append(prominences[i]*prominences[i+n])

    return bins, distances, weights

vert_bin, vert_dist, vert_weight = calculate_distances_to_n_nearest_neighbors_of_1D_array(vert_max[0], vert_max[1]["prominences"],15)
#vert_mean_bins = calculate_distances_to_n_nearest_neighbors_of_1D_array(vert_max, 10)[1]
print(vert_dist)

# plot histogram of vertical distances
# new plot window: 
plt.figure()
import matplotlib.pyplot as plt
plt.hist(vert_dist, bins = 500)

# if perfectly calibrated, the delta m/z should be peaked in a distance of 1.007276466812
# do an pre-guessed calibration, e.g. by TOFWerk 
a=0.0
b=2

vert_dist_pre_calib = [a+x**b for x in vert_dist]

plt.figure()
plt.hist(vert_dist_pre_calib, bins = 300)

plt.figure()
# plt.plot(vert_bin, vert_dist, "x", color = "red", alpha = vert_weight/np.max(vert_weight))
# plot vert_bin vs vert_dist with color according to vert_weight
plt.scatter(vert_bin, vert_dist, c = vert_weight, cmap = "turbo", alpha = 0.5)

# %%
#masscalibration_p1 = 2829.4940436162415
#masscalibration_p2 = -4042.695405928985
#masscalibration_p3 = 0.49993765722959954

vert_bin, vert_dist, vert_weight = calculate_distances_to_n_nearest_neighbors_of_1D_array((vert_max[0]**1/masscalibration_p3+masscalibration_p2)/masscalibration_p1, vert_max[1]["prominences"],3)

# plot histogram of vertical distances in mass space
# new plot window: 
plt.figure()
import matplotlib.pyplot as plt
plt.hist(vert_dist, bins = 5000)

plt.figure()
# plt.plot(vert_bin, vert_dist, "x", color = "red", alpha = vert_weight/np.max(vert_weight))
# plot vert_bin vs vert_dist with color according to vert_weight
plt.scatter(vert_bin, vert_dist, c = vert_weight, cmap = "turbo", alpha = 0.5)

# retention time 

#hori_dist = calculate_distances_to_n_nearest_neighbors_of_1D_array(hori_max[0], 10)
#print(hori_dist)

## plot histogram of vertical distances
#import matplotlib.pyplot as plt
#plt.figure()
#plt.hist(hori_dist, bins = 1000)


# %%

def calculate_distances_of_entries_of_1D_array(array):
    distances = []
    for i in range(len(array)-1):
        distances.append(array[i+1]-array[i])
    return distances

hori_dist = calculate_distances_of_entries_of_1D_array(hori_max[0])

# %%
def plot_histogram_of_vertmax(hori_dist):
    import matplotlib.pyplot as plt
    plt.hist(hori_dist)
    plt.show()

plot_histogram_of_vertmax(vert_dist)



# %%
# CORRECT FOR BACKGROUND










# %%
# EVALUATE VALUES OF GRIDD CROSSINGS
def evaluate_values_at_gridded_peaks(filtered_im_as_np, vert_max, hori_max):
    # define grid
    grid = np.zeros_like(filtered_im_as_np)
    for i in vert_max[0]:
        for j in hori_max[0]:
            grid[j,i] = 1
    # find maxima in grid
    grid_max = find_max_2D(grid)
    # find maxima in image
    image_max = find_max_2D(filtered_im_as_np)
    # find maxima in image and grid
    grid_image_max = np.logical_and(grid_max, image_max)
    # find values at maxima in image and grid
    grid_image_max_values = filtered_im_as_np[grid_image_max]
    return grid_image_max_values

vals_on_grid = evaluate_values_at_gridded_peaks(filtered_im_as_np, vert_max, hori_max) 

print(vals_on_grid)


# FIND MAXIMA IN 2D IMAGE
def define_gridded_peaks(filtered_im_as_np, vert_max, hori_max):
    # define grid
    grid = np.zeros_like(filtered_im_as_np)
    for i in vert_max[0]:
        for j in hori_max[0]:
            grid[j,i] = 1
    # find maxima in grid
    grid_max = find_max_2D(grid)
    # find maxima in image
    image_max = find_max_2D(filtered_im_as_np)
    # find maxima in image and grid
    grid_image_max = np.logical_and(grid_max, image_max)
    return grid_image_max

# %%
# CONVERT MAXIMA TO LIST OF X,Y,Z including only values larger than zero
threshold = 0.5
def convert_matrix_elements_larger_than_zero_to_list_x_y_z(matrix, filtered_im_as_np):
    # convert matrix elements larger than zero to coordinate
    matrix = matrix.astype(int)
    matrix = np.where(matrix > threshold, matrix, 0)
    matrix = np.array(matrix).T
    # convert to list of x,y,z
    list_x_y_z = []
    for i in range(len(matrix)):
        list_x_y_z.append([matrix[i][0], matrix[i][1], filtered_im_as_np[matrix[i][0], matrix[i][1]]])
    return np.array(list_x_y_z)


grid_image_max = define_gridded_peaks(filtered_im_as_np, vert_max, hori_max)
coordinates_list = convert_matrix_elements_larger_than_zero_to_list_x_y_z(grid_image_max, filtered_im_as_np)

def color_plot_x_y_z_list(x_y_z_list):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #ax.scatter(x_y_z_list[:,0], x_y_z_list[:,1], x_y_z_list[:,2])
    ax.stem(x_y_z_list[:,0], x_y_z_list[:,1],  x_y_z_list[:,2])
    plt.show()

color_plot_x_y_z_list(coordinates_list)


def plot_grid_image_max(grid_image_max):
    import matplotlib.pyplot as plt
    plt.imshow(grid_image_max)
    plt.show()

# %%

def plot_stem_3D(x_y_z_list):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # add horizontal lines
    for i in hori_max[0]:
        ax.plot([i,i], [0, np.max(vert_max[0])], [0,0], 'k', linewidth=0.5)
    # add vertical lines
    for i in vert_max[0]:
        ax.plot([0, np.max(hori_max[0])], [i,i], [0,0], 'k',linewidth=0.5)

    ax.stem(x_y_z_list[:,0], x_y_z_list[:,1], x_y_z_list[:,2], linefmt='C0-', markerfmt='C0o', basefmt=" ")
    plt.show()

plot_stem_3D(coordinates_list)
 # %% PLOT NDARRAY AS IMAGE WITH RECOGNIZED PEAKS

def plot_ndarray_as_image_with_colorbar_with_coordinatelist_as_scatterplot(ndarray, coordinates_list):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # add horizontal lines
    for i in hori_max[0]:
        ax.plot( [0, np.max(vert_max[0])],[i,i],'k', linewidth=0.5)
    # add vertical lines
    for i in vert_max[0]:
        ax.plot( [i,i],[0, np.max(hori_max[0])], 'k',linewidth=0.5)
    # plot image
    ax.imshow(ndarray, cmap='gray', origin='lower')
    # plot scatterplot
    ax.scatter(coordinates_list[:,1], coordinates_list[:,0], coordinates_list[:,2], edgecolors='r', facecolors = 'none')
    plt.show()

plot_ndarray_as_image_with_colorbar_with_coordinatelist_as_scatterplot(filtered_im_as_np, coordinates_list)
# %%

def plot_image_as_image(im_as_np):
    import matplotlib.pyplot as plt
    plt.imshow(im_as_np)
    plt.show()

plot_image_as_image(filtered_im_as_np)
# %%



def plot_image_including_grid_and_coordinate_list_z_as_color(filtered_im_as_np, coordinates_list, vert_max, hori_max):
    import matplotlib.pyplot as plt
    from matplotlib import cm   
    #from mpl_toolkits.mplot3d.axes3d import plot_surface
    from mpl_toolkits import mplot3d
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # add horizontal lines
    for i in hori_max[0]:
        ax.plot([i,i], [0, np.max(vert_max[0])], [0,0], 'k', linewidth=0.5)
    # add vertical lines
    for i in vert_max[0]:
        ax.plot([0, np.max(hori_max[0])], [i,i], [0,0], 'k',linewidth=0.5)
    # add image
    lenx =len(filtered_im_as_np[0])
    leny =len(filtered_im_as_np[1])

    # get x, y, z coordinates from filtered_im_as_np
    X, Y = np.meshgrid(range(lenx), range(leny))
    Z = filtered_im_as_np



    ax.plot_surface(X, Y, Z, cmap='gray')
    # add coordinates
    ax.scatter(coordinates_list[:,0], coordinates_list[:,1], coordinates_list[:,2], c='r')
    plt.show()

plot_image_including_grid_and_coordinate_list_z_as_color(filtered_im_as_np, coordinates_list, vert_max, hori_max)

#%%


#print dimensions of coordinates_list
print("dimensions of coordinates_list")
print(coordinates_list.shape)
print(coordinates_list)





print(grid_image_max.shape)
print(coordinates)
print(coordinates.shape)



def plot_grid_image_max(grid_image_max):
    import matplotlib.pyplot as plt
    plt.imshow(grid_image_max)
    plt.show()

plot_grid_image_max(grid_image_max)





# %%
