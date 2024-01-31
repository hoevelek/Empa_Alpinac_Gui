#%%
#import sys
#sys.path.append('c:\\Users\\kaho\\polybox\\ALPINAC\\Alpinac_packing\\Alpinac_packing\\alpinac_sups')
import time
t_start = time.time()
from pathlib import Path

from alpinac_sups.read_h5 import read_desc, read_ionisation_modes, read_mass_axis, H5_Reader, read_time_axis
from matplotlib.lines import segment_hits



#
#file = Path(r"C:\Users\coli\Documents\TW_EI_CI_data\DataFile_2021.09.26-16h56m06s.h5")


file = Path(r"C:\Users\kaho\TOFWerk_data\shared_data\ExtractTofSpectra_example\lib\DataFile_2021.09.26-16h56m06s.h5")
file = Path(r"G:\503_Themen\Klima\DNAtmos\test_data221115\DataFile_2021.09.26-16h56m06s.h5")

ei_seg = (1,5)
ci_seg = (6,9)


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

#ionisation_dict = {
#    # Last value is not included
#    "EI": ei_seg,
#    "CI": ci_seg,
#}




# %% Reader class

reader = H5_Reader(

    file,
    ionisation_dict

)

# %% This will read only the EI data

data = reader.read_full_ionization_data('EI')

# %% This will return the whole data over all segements

data_all = reader.read_full_tof_data()

data[0,1]

type(data)



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


import numpy as np
from PIL import Image, ImageOps, ImageFilter
from scipy.ndimage import maximum_filter, gaussian_filter


    # FILTER TO DETECT VERTICL PLUS HORIZONTAL STEPS AND RETURN BINARY=1 IF DETECTED, 0 ELSE
    #%%



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


# EXTRACT PEAKS

#convert to 8bit image
im_raw = Image.fromarray(data).convert("L")

 # SMOOTH IN 2D
im_autocontrast_filter = im_raw.filter(ImageFilter.SMOOTH)#.filter(ImageFilter.FIND_EDGES)

im_autocontrast_filter.save(Path(r"C:\Users\kaho\data_Empa\210623.1218.std.9.2Dpic_autocontrast_smoothed.jpeg"))

# CONVERT BACK TO NP ARRAY (REMOVE OTHER ARRAYS FROM MEM IF READY!)
filtered_im_as_np = np.asarray(im_autocontrast_filter)

# GENERATE RT AND MASS SUM SPECTRA OF INT8 IMAGE DATA - SMOOTHED

# find stripes, all is binary
vertical_lines = np.sum(filtered_im_as_np[0:int(dim_pic[0]/1)], axis = 0)
horizontal_lines = np.sum(filtered_im_as_np[0:int(dim_pic[0]/1)], axis = 1)
max_vert = np.max(vertical_lines)

from scipy.signal import find_peaks
vert_max = find_peaks(vertical_lines, prominence = (max_vert*0.01, None), distance =7)
hori_max = find_peaks(horizontal_lines)

