from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
from scipy import special
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn.datasets import make_blobs
from alpinac_sups.read_h5 import read_desc, read_ionisation_modes, read_mass_axis, H5_Reader, read_time_axis
# find_peaks_cwt
from scipy.signal import find_peaks_cwt, find_peaks
from scipy.optimize import curve_fit
from scipy.optimize import optimize as optim
from scipy.optimize import minimize
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
def sum_over_dim(data, dim, sum_dim, sum_range):
    # data: data array
    # dim: dimension to sum over
    # sum_range: range of dim to sum over
    # returns: sum of data over dim in range sum_range
    if dim == 0:
        return np.sum(data[sum_range[0]:sum_range[1],:], axis=sum_dim)
    elif dim == 1:
        return np.sum(data[:,sum_range[0]:sum_range[1]], axis=sum_dim)
    else:
        print('dim must be 0 or 1')






# Funktion to plot one dimension of the data by summing over selected part of the other dimension
def plot_1D(data, dim, sum_dim, sum_range):
    # data: data array
    # dim: dimension to plot
    # sum_dim: dimension to sum over
    # sum_range: range of sum_dim to sum over
    # returns: plot of 1D data
    fig, ax = plt.subplots()
    plot_data = sum_over_dim(data, dim, sum_dim, sum_range)
    plt.plot(plot_data)
    plt.show()

#retention time plot
plot_1D(data=data, dim=0, sum_dim=1, sum_range=(0, data.shape[1]))

#mass plot
plot_1D(data=data, dim=1, sum_dim=0, sum_range=(0, data.shape[0]))

# get the average noise. Because most of the spectrum has no signal, the median of rt column-wise stadard deviations is a good estimate for the noise

# get the standard deviation of each column:
std = np.std(data, axis=0)
# get the median of the standard deviations:
noise = np.median(std)
print('noise: ', noise)

# it is likely that everything above 2*noise is not possible to fit. So we can use this as a threshold for the peakfinder.
# We can also use this to filter out the noise from the data. This is done by replacing all values above 2*noise with 0.
mask = data < 3*noise
data_filtered = np.where(mask, 0, data)
print(data_filtered.shape)
print(data.shape)

plot_1D(data=data, dim=0, sum_dim=1, sum_range=(0, data.shape[1]))
plot_1D(data=data_filtered, dim=0, sum_dim=1, sum_range=(0, data_filtered.shape[1]))




#data_filtered = np.where(data > 2*noise, 0, data)
# plot filtered data sum over rt togehter with the original data
fig, ax = plt.subplots()
plot_data = sum_over_dim(data=data, dim=0, sum_dim=1, sum_range=(0, data.shape[1]))
plot_data_filtered = sum_over_dim(data=data_filtered, dim=0, sum_dim=1, sum_range=(0, data_filtered.shape[1]))
plt.plot(plot_data)
plt.plot(plot_data_filtered)
plt.show()


#data_filtered = np.where(data > 2*noise, 0, data)
# plot filtered data sum over rt togehter with the original data
fig, ax = plt.subplots()
plot_data = sum_over_dim(data=data, dim=1, sum_dim=0, sum_range=(0, data.shape[0]))
plot_data_filtered = sum_over_dim(data=data_filtered, dim=1, sum_dim=0, sum_range=(0, data_filtered.shape[0]))
plt.plot(plot_data)
plt.plot(plot_data_filtered)
plt.show()

# Try to seperate co-eluting data

# try to find contributing masses:
fig, ax = plt.subplots()
mass_spectra_for_certain_rt = sum_over_dim(data_filtered[6050:6350,:], dim=1, sum_dim=0, sum_range=(0,data_filtered.shape[0]))
plt.plot(mass_spectra_for_certain_rt)
mass_spectra_for_certain_rt = sum_over_dim(data_filtered[6050:6200,:], dim=1, sum_dim=0, sum_range=(0,data_filtered.shape[0]))
plt.plot(mass_spectra_for_certain_rt)
mass_spectra_for_certain_rt = sum_over_dim(data_filtered[6200:6350,:], dim=1, sum_dim=0, sum_range=(0,data_filtered.shape[0]))
plt.plot(mass_spectra_for_certain_rt)
plt.show()




fig, ax = plt.subplots()
p1 = sum_over_dim(data=data_filtered[:,14110:14130], dim=0, sum_dim=1, sum_range=(0, data_filtered.shape[1]))
p2 = sum_over_dim(data=data_filtered[:,13885:13889], dim=0, sum_dim=1, sum_range=(0, data_filtered.shape[1]))
plt.plot(p1)
plt.plot(p2)
plt.show

# plot relative contribution of each mass to the total signal
fig, ax = plt.subplots()
plt.plot(p2/p1)
plt.show()

# sum of 3 gaussian peaks
def gauss(x, a, b, c, d, e, f, g, h, i):
    return a*np.exp(-((x-b)/c)**2) + d*np.exp(-((x-e)/f)**2) + g*np.exp(-((x-h)/i)**2)



sigma = 1.0
gamma = 2.0
fraction = 0.5

from scipy.special import voigt_profile 

# sum of pseudo voigt peaks
def pseudo_voigt(x, a, b, c, d, e, f, g, h, i, k, l, m):
    return a*voigt_profile(x-b, c, d) + e*voigt_profile(x-f, g, h) + i*voigt_profile(x-k, l, m)


# plot the pseudo voigt peaks

fig, ax = plt.subplots()
x = np.linspace(0, 20, 1000)
plt.plot(x, pseudo_voigt(x, 1, 10, 0.3, 1.1, 1, 12, 0.3, 1.1, 1, 14, 0.3, 1.1))
plt.plot(x, pseudo_voigt(x, 1, 10, 0.4, 1.1, 1, 12, 0.4, 1.1, 1, 14, 0.4, 1.1))

plt.show()

fig, ax = plt.subplots()
# same heigth, same width, different position
same = plt.plot( x,  pseudo_voigt(x, 1, 10, 0.3, 1.1, 1, 12, 0.3, 1.1, 1, 14, 0.3, 1.1)/pseudo_voigt(x, 1, 10, 0.3, 1.1, 1, 12, 0.3, 1.1, 1, 14, 0.3, 1.1))
first_low = plt.plot( x,  pseudo_voigt(x, 0.2, 10, 0.3, 1.1, 1, 12, 0.3, 1.1, 1, 14, 0.3, 1.1)/pseudo_voigt(x, 1, 10, 0.3, 1.1, 1, 12, 0.3, 1.1, 1, 14, 0.3, 1.1))
middle_low = plt.plot( x,  pseudo_voigt(x, 1, 10, 0.3, 1.1, 0.2, 12, 0.3, 1.1, 1, 14, 0.3, 1.1)/pseudo_voigt(x, 1, 10, 0.3, 1.1, 1, 12, 0.3, 1.1, 1, 14, 0.3, 1.1))
last_low = plt.plot( x,  pseudo_voigt(x, 1, 10, 0.3, 1.1, 1, 12, 0.3, 1.1, 0.2, 14, 0.3, 1.1)/pseudo_voigt(x, 1, 10, 0.3, 1.1, 1, 12, 0.3, 1.1, 1, 14, 0.3, 1.1))
first_high = plt.plot( x,  pseudo_voigt(x, 1, 10, 0.3, 1.1, 0.2, 12, 0.3, 1.1, 0.2, 14, 0.3, 1.1)/pseudo_voigt(x, 1, 10, 0.3, 1.1, 1, 12, 0.3, 1.1, 1, 14, 0.3, 1.1))
middle_high = plt.plot( x,  pseudo_voigt(x, 0.2, 10, 0.3, 1.1, 1, 12, 0.3, 1.1, 0.2, 14, 0.3, 1.1)/pseudo_voigt(x, 1, 10, 0.3, 1.1, 1, 12, 0.3, 1.1, 1, 14, 0.3, 1.1))
last_high = plt.plot( x,  pseudo_voigt(x, 0.2, 10, 0.3, 1.1, 0.2, 12, 0.3, 1.1, 1, 14, 0.3, 1.1)/pseudo_voigt(x, 1, 10, 0.3, 1.1, 1, 12, 0.3, 1.1, 1, 14, 0.3, 1.1))
# add legend
plt.legend((same[0], first_low[0], middle_low[0], last_low[0], first_high[0], middle_high[0], last_high[0]), ('same', 'first_low', 'middle_low', 'last_low', 'first_high', 'middle_high', 'last_high'))
#last_low = plt.plot( x,  pseudo_voigt(x, 1, 10, 0.9, 1.1, 1, 12, 0.9, 1.1, 0.2, 14, 0.9, 1.1)/pseudo_voigt(x, 1, 10, 0.9, 1.1, 1, 12, 0.9, 1.1, 1, 14, 0.9, 1.1))
# plot vertical lines at 10,12,14
plt.axvline(x=10, color='k', linestyle='--')
plt.axvline(x=12, color='k', linestyle='--')
plt.axvline(x=14, color='k', linestyle='--')
plt.show()

# for 3 co-eluting substances there are 6 classes of peak ratios



#
# determine all maxima in the mass spectra

# get the mass spectra for certain rt

mass_spectra_for_certain_rt = sum_over_dim(data_filtered[6050:6350,:], dim=1, sum_dim=0, sum_range=(0,data_filtered.shape[0]))

# get the maxima using find_peaks



maxima = find_peaks_cwt(mass_spectra_for_certain_rt, widths = 3, min_snr=3)
# get maxima positions
maxima = maxima


# calculate the differenece between each 6 neighboured maxima positions and plot it as a histogram
maxima = find_peaks_cwt(mass_spectra_for_certain_rt, widths = 3, min_snr=noise)
# get the difference between each 6 neighboured maxima
diff = np.diff(maxima**2, n=1)
# filter out the differences that are smaller than 100**2
diff = diff[diff>150**2]
# plot the histogram
fig, ax = plt.subplots()
plt.hist(diff, bins=1000)
plt.show()



#  AUTO MASS CALIBRATION

#convert histogram to xy data
#fit pseudo voigt function to the histogram
# get appropriate bins
bins = np.linspace(0, max(diff), 1000)
# get the histogram bins and counts


hist, bins = np.histogram(diff, bins=bins)
# get the bin centers
bin_centers = (bins[:-1] + bins[1:])/2

# get appropriate bins
fig, ax = plt.subplots()
plt.plot(bin_centers, hist)
plt.show()




def voigt_profile_fit(x, A, mu, sigma, gamma):
    return A*voigt_profile(x-mu, sigma, gamma)  


def optim_func(params, x):
    #a=x0[0]
    print(x)
    scale = params['scale']*105
    power = params['power']*0.5
    shift = params['shift']*400

    print(params['shift'])
    print(params['power'])
    print(params['scale'])
    diffr = 11 - np.diff(10 + (scale*(x-shift))**power, n=1)
    # bins = np.linspace(0, np.max(diff), 1000)
    # # get the histogram bins and counts
    # hist, bins = np.histogram(diff, bins=bins)
    # # get the bin centers
    # bin_centers = (bins[:-1] + bins[1:])/2
    # hist_filt = savgol_filter(hist, 51, 3)
    # #get maximum value
    # max_val = np.max(hist_filt)
    # # get first bin_centers value where hist_filt is larger than 0.5
    # first = bin_centers[np.where(hist_filt>0.5*max_val)[0][0]]
    # # get last bin_centers value where hist_filt is larger than 0.5
    # last = bin_centers[np.where(hist_filt>0.5*max_val)[0][-1]]
    # width = last - first
    # # calculate the FWHM of the bin_centers and hist_filter curve
    # # FWHM numerical
    # # get the largest maximum of the histogram
    # # maxima = find_peaks(hist_filt, width = 2e3)[0]
    # # get the maxima positions and the maxima values
    # # find the maxima positions
    # popt, pcov = curve_fit(voigt_profile_fit, bin_centers, hist, p0=[6e6,1.8e5,3.66000131e+01,9.38346658e+04])
    # A, mu, sigma, gamma = popt
    # width = 0.5343*gamma + np.sqrt(0.2169*gamma**2 + sigma**2)
    # #height = A
    # print(width)
    # print(np.std(diff-1)/np.mean(diff))
    # print(np.std(diff)/np.mean(diff))
    # print(mu)
    print(np.mean(diffr))
    print(np.mean(diffr)-10)
    ret = np.mean(diffr)
    print(diffr[0:5])
    return np.mean(diffr)

from lmfit import Parameters, Minimizer
guess_params = Parameters()
guess_params.add('scale', value = float(2800/2826), min = 0, max = 5)
guess_params.add('power', value = float(0.5/0.5005),min = 0.5*1.5, max=0.5*2.5)
guess_params.add('shift', value = float(-4015/4000), min = -5, max = 5)
#guess_params.add('p5', value = 0, vary = False)

minner = Minimizer(optim_func, guess_params, fcn_args=[x]) #, method='leastsq'
tof_to_mass_optimised = minner.minimize(method = 'nelder', params = guess_params)

scale = guess_params['scale']*105
power = guess_params['power']*0.5
shift = guess_params['shift']*400
x = maxima
np.mean(1 - np.diff((scale*abs(x-shift))**power, n=1))




fig, ax = plt.subplots()
params = guess_params
np.mean(np.diff(3.958763280972075/2800*(x-0.29262057364255867*4000)**(0.9997631991683984*0.5), n=1)-[1]*len(diff))
diff = np.diff(params['scale']/100*(x-params['shift'])**params['power'], n=1)
bins = np.linspace(0, np.max(diff), 1000)
# get the histogram bins and counts
hist, bins = np.histogram(diff, bins=bins) 
#plot histogram
plt.hist(diff, bins)

plt.hist()


x = maxima

# minimize the width of the pseudo-voigt profile by optimizing a,b,c
res = minimize(optim_func, x0=[2.0,-20], args=(x), method='Nelder-Mead', tol=1e-8)
print(res.x)

#plot the histogram and the fit
fig, ax = plt.subplots()

#rolling avg of hist to smooth it
import numpy as np
from scipy.signal import savgol_filter
savgol_filter(hist, 51, 3)


plt.plot(bin_centers, hist)
plt.plot(bin_centers, savgol_filter(hist,51,3))
diff_fit =  res.x[0]*(x**res.x[1]) + res.x[2]
bins_fit = np.linspace(0, max(diff_fit), 1000)
# get the histogram bins and counts
hist_fit_plot, bin_fit_plot = np.histogram(diff_fit, bins=bins_fit)
# get the bin centers
bin_centers_fit_plot = (bin_fit_plot[:-1] + bin_fit_plot[1:])/2
plt.plot(bin_centers_fit_plot, hist_fit_plot)
plt.show()





# fit the histogram
# plot the histogram and the fit
fig, ax = plt.subplots()
plt.plot(bin_centers, hist)
plt.plot(bin_centers, voigt_profile_fit(bin_centers, *popt))
plt.show()
np.array([6.16446974e+06, 1.72517621e+05, 3.66000131e+01, 9.38346658e+04])
width = 0.5343*gamma + np.sqrt(0.2169*gamma**2 + sigma**2)
heigth = A
# do an optimization to result in the pseudo-voigt profile with the heighest amplitude by varying a,b of diff = a*(maxima**b) + c

# define the function to be optimized
def func(x, a, b, c):
    return a*(x**b) + c

# get the difference between each 6 neighboured maxima





# plot mass spectra and maxima
fig, ax = plt.subplots()
plt.plot(mass_spectra_for_certain_rt)
plt.plot(maxima, mass_spectra_for_certain_rt[maxima], "x")
plt.show()





# plot the gaussian peaks
fig, ax = plt.subplots()
x = np.linspace(0, 15, 1000)
plt.plot(x, gauss(x, 1, 10, 0.3, 1, 11, 0.3, 1, 12, 0.3))
plt.show()

#fig, ax = plt.subplots()
plt.plot( x,  gauss(x, 0.3, 10, 0.3, 1, 11, 0.3, 1, 12, 0.3)/gauss(x, 1, 10, 0.3, 0.3, 11, 0.3, 0.7, 12, 0.3))
plt.show()










# plot the standard deviation histogram
# get appropriate bins
bins = np.linspace(0, 20*noise, 1000)
# plot
fig, ax = plt.subplots()
hist=np.histogram(std, bins=bins)
plt.plot(hist[1][:-1], hist[0])
plt.show()







# calc


# seperate co-eluting peaks
# use peakfinder to find peaks in retention time plot
from scipy.signal import find_peaks as fp
peaks, _ = fp(plot_data, height=1000, distance=100)


