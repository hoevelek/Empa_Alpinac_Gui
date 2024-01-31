import matchms
import numpy as np

from matchms import Spectrum
from matchms.filtering import normalize_intensities
# import formula package
from molmass import Formula as chem_formula



# generate random spectrum with mz and intensities
np.random.seed(42)
number_of_peaks = 10
mz = np.random.uniform(low=0, high=300, size=number_of_peaks)
mz = np.array([20, 22, 24, 26, 28, 30, 32, 34, 36, 38]).astype(float)
mz_order = np.argsort(mz)
mz = mz[mz_order]
intensities = np.random.uniform(low=0.01, high=1, size=number_of_peaks)
intensities = intensities[mz_order]
spectrum_exp = Spectrum(mz=mz, intensities=intensities)
# Normalize intensities
spectrum_exp = normalize_intensities(spectrum_exp)


# load spectrum from jdx file
path = r"C:\Users\kaho\polybox\DNAtmos_code\tools_for_Empa\594-18-3-Mass.jdx" #CBr2Cl2
path = r"C:\Users\kaho\polybox\DNAtmos_code\tools_for_Empa\630-25-1-Mass.jdx" #C2Br2Cl4

# read NIST data
from alpinac.io_tools import get_data_from_jdx_file
from jcamp import jcamp_read #to read NIST spectra

with open(path, 'r') as f:
    dict_NIST = jcamp_read(f)

# get mz and intensities
mz = dict_NIST['x']
intensities = dict_NIST['y']
# create spectrum
spectrum_exp = Spectrum(mz=mz, intensities=intensities)
# Normalize intensities
spectrum_exp = normalize_intensities(spectrum_exp)

# read in table from C:\Users\kaho\polybox\DNAtmos_code\tools_for_Empa
path_rel_spectra = r"C:\Users\kaho\polybox\DNAtmos_code\tools_for_Empa\ClBr.csv"
# read in table from C:\Users\kaho\polybox\DNAtmos_code\tools_for_Empa as pandas
import pandas as pd

df = pd.read_csv(path_rel_spectra, sep=",", header=0)
# generate first relative spectrum
# mz are from second column on the parts after "+" of the strings in the first line
mz = df.columns[1:].str.split("+", expand=True)
# get astype function from numpy
from numpy import array

# get second element of each list element of mz
mz = np.array([mz_elem[1] for mz_elem in mz]).astype(float)
# get formula names of first columns
formula_names = df.iloc[:, 0]
# for each formula name get the relative intensities, and use mz to create matchms.Spectrum for each formula and name it by the formula name
spectra_rel_initial = []
for i, formula_name in enumerate(formula_names):
    # get intensities
    intensities = np.array(df.iloc[i, 1:]).astype(float)
    # get indici of intensities that are a float > 0
    valid_indices = np.where(intensities > 0)
    # get mz and intensities of the relative spectrum
    mz_rel_i = np.array(mz)[valid_indices]
    intensities_rel_i = np.array(intensities)[valid_indices]
    # create spectrum
    spectra_rel_initial.append(
        Spectrum(mz=mz_rel_i, intensities=intensities_rel_i, metadata={"formula": formula_names[i]})
    )


# function spectrum shifter: shifts the mz values of a spectrum by a given value
def spectrum_shifter(spectrum, shift):
    # create new spectrum with shifted mz values and same intensities
    mz_new = spectrum.peaks.mz + shift
    spectrum_out = Spectrum(mz=mz_new, intensities=spectrum.peaks.intensities)
    # add shift to metadata
    spectrum_out.set("shift", shift)
    # add formula to metadata
    spectrum_out.set("formula", spectrum.get("formula"))
    return spectrum_out


def get_part_of_exp_spectrum(mz_exp, mz_rel):
    # create a spectrum of the mz values of the relative spectrum and intensities of the experimental spectrum at those mz values. If mz value is not in mz_exp, set intensity to 0
    spectrum_exp_part = Spectrum(
        mz=mz_rel,
        intensities=np.array(
            [
                spectrum_exp.peaks.intensities[np.where(mz_exp == mz_rel_i)]
                if mz_rel_i in mz_exp
                else 0
                for mz_rel_i in mz_rel
            ]
        ).astype(float).flatten(),
    )
    return spectrum_exp_part



# test spectrum shifter
#spectrum_shifted = spectrum_shifter(spectrum_exp, 10)
# get all masses of the experimental spectrum
mz_exp = spectrum_exp.peaks.mz

# creates spectra with shifts from 1:300
spectra_rel_shifted = []
for spectrum in spectra_rel_initial:
    spectra_rel_shifted_part = [spectrum_shifter(spectrum, i) for i in mz_exp]
    spectra_rel_shifted.extend(spectra_rel_shifted_part)

# get 5 highest matching spectra and printout their formula names and shifts
from matchms.similarity import CosineGreedy

# initialize similarity measure
similarity_measure = CosineGreedy(tolerance=0.1, mz_power=5, intensity_power=1.5)

# calculate similarity between shifted spectrum and all relative spectra
scores = [
    similarity_measure.pair(reference=spectrum_rel_shifted, query=spectrum_exp)
    for spectrum_rel_shifted in spectra_rel_shifted
]

scores = []
for spectrum_rel_shifted in spectra_rel_shifted:
    # get part of experimental spectrum that contains the masses where the relative spectrum has intensities
    spectrum_exp_part = get_part_of_exp_spectrum(spectrum_exp.peaks.mz, spectrum_rel_shifted.peaks.mz)
    if sum(spectrum_exp_part.peaks.intensities) < 1E-18:
        score = 0
    else:
    # calculate similarity between shifted spectrum and all relative spectra
        score = similarity_measure.pair(query=spectrum_rel_shifted, reference=spectrum_exp_part)
    scores.append(score)
    


#calculate a list with start m/z, end m/z and difference of m/z for all m/z values of the experimental spectrum
mz_exp_diff = np.diff(mz_exp)
# calculate diff from all to all other mz values
mz_exp_diff = np.array([mz_exp_i - mz_exp for mz_exp_i in mz_exp])#.flatten()
# get all intensities of the experimental spectrum
weight_exp_diff = np.array([spectrum_exp.peaks.intensities[np.where(mz_exp == mz_exp_i)] * spectrum_exp.peaks.intensities for mz_exp_i in mz_exp])#.flatten()
# get unique mz_exp_diff values > 0 
#mz_exp_diff = mz_exp_diff[mz_exp_diff > 0]
# get unique mz_exp_diff values > 0
mz_uniques = np.unique(mz_exp_diff)
# do pd frame with mz_exp_diff, weight_exp_diff and occurrences
dict_weighted={}
for mz_unique in mz_uniques:
    if mz_unique > 0:
        mask_mz_unique = mz_exp_diff == mz_unique
        occurrences = np.sum(mask_mz_unique)
        # get number of occurrences
        weighted_occurences = np.sum(weight_exp_diff[mask_mz_unique]*mz_exp_diff[mask_mz_unique])
        # get weighted occurrences
        dict_weighted[mz_unique] = {"occurrences": occurrences, "weighted_occurences": weighted_occurences}

# convert dict to pd frame
df_weighted = pd.DataFrame.from_dict(dict_weighted, orient="index")

# add barplot of chlorinated differences:
# get all mass diffs whcih are combinations of n times 35 and n2 time 37 and < 300
combis_chlorine = [35,36,37,38,49,50, 51, 52, 69, 83, 84, 91, 93] # 84: LM
for i in range(0, 10):
    for j in range(0, 10):
        comb = 35*i + 37*j
        if 0 < comb < 300:
            combis_chlorine.append(comb)

combis_bromine = [79, 80, 81, 82, 93, 95, 135, 137]
for i in range(0, 10):
    for j in range(0, 10):
        comb = 79*i + 81*j
        if 0 < comb < 300:
            combis_bromine.append(comb)

combis_fluorine = [19,20,33]
for i in range(0, 20):
        comb = 19*i
        if 0 < comb < 300:
            combis_fluorine.append(comb)


#combis_carbons = []
#for i in range(0, 20):
#        comb = 12*i 
#        if 0 < comb < 300:
#            combis_carbons.append(comb)

combis_hydrocarbons = [26, 27, 28, 29, 41, 42, 43, 50, 51, 52, 56, 57, 65, 71, 76, 77, 78, 85, 91, 92, 103, 104, 105, 115]

combis_iodine = [127, 128, 142, 156]

pollutions = [44]




#combis_hydrocarbons.append(1)
# combis_hydrocarbons = []
# for i in range(0, 20):
#     for j in range(0, 20):
#         comb = 1*i + 2*j
#         if comb < 300:
#             combis_hydrocarbons.append(comb)
        
# get intensities of the experimental combis from the histogram
#intensities_chlorine = np.array([df_weighted["occurrences"][df_weighted.index == combi] for combi in combis_chlorine]).flatten()



# plot histograms of occurrences and weighted occurrences
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1,1)
plt.bar(df_weighted.index, df_weighted["occurrences"], width = 0.5, color='grey', alpha=0.5)
plt.ylabel('Number of peaks')
plt.xlabel('m/z difference')
plt.show()

import matplotlib.pyplot as plt
fig, ax = plt.subplots(1,1)
plt.bar(df_weighted.index, df_weighted["weighted_occurences"], width = 0.5, color='grey', alpha=0.5)
plt.ylabel('Number of peaks')
plt.xlabel('m/z difference')
plt.show()

# add chlorinated differences using dict_weighted and plot them
#intensities_carbons = [dict_weighted[ind_combi]["weighted_occurences"] if ind_combi in dict_weighted.keys() else 0 for ind_combi in combis_carbons]
#plt.bar(combis_carbons, intensities_carbons, width = 0.5, color='black', alpha=0.5)

intensities_hydrocarbons = [dict_weighted[ind_combi]["weighted_occurences"] if ind_combi in dict_weighted.keys() else 0 for ind_combi in combis_hydrocarbons]
plt.bar(combis_hydrocarbons, intensities_hydrocarbons, width = 0.5, color='black', alpha=0.5)

intensities_fluorine = [dict_weighted[ind_combi]["weighted_occurences"] if ind_combi in dict_weighted.keys() else 0 for ind_combi in combis_fluorine]
plt.bar(combis_fluorine, intensities_fluorine, width = 0.5, color='green', alpha=0.5)

intensities_chlorine = [dict_weighted[ind_combi]["weighted_occurences"] if ind_combi in dict_weighted.keys() else 0 for ind_combi in combis_chlorine]
plt.bar(combis_chlorine, intensities_chlorine, width = 0.5, color='purple', alpha=0.5)

intensities_bromine = [dict_weighted[ind_combi]["weighted_occurences"] if ind_combi in dict_weighted.keys() else 0 for ind_combi in combis_bromine]
plt.bar(combis_bromine, intensities_bromine, width = 0.5, color='red', alpha=0.5)

intensities_iodine = [dict_weighted[ind_combi]["weighted_occurences"] if ind_combi in dict_weighted.keys() else 0 for ind_combi in combis_iodine]
plt.bar(combis_iodine, intensities_iodine, width = 0.5, color='blue', alpha=0.5)














# # do a histogram of mz_exp_diff and correlating weightfactors (multiplied for all contributing)
# hist, bin_edges = np.histogram(mz_exp_diff, bins = int(max(mz_exp_diff)-min(mz_exp_diff)))
# # get weights by masking the intensities of the experimental spectrum with the mz_exp_diff
# mask = np.isin(mz_exp_diff, bin_edges)
# # intensity matrix
# intensity_matrix = np.array([spectrum_exp.peaks.intensities[np.where(mz_exp == mz_exp_i)] for mz_exp_i in mz_exp])
# #mz_exp_diff = mz_exp_diff[mz_exp_diff > 9]
# # get histogram of mz_exp_diff
# hist, bin_edges = np.histogram(mz_exp_diff, bins = int(max(mz_exp_diff)-min(abs(mz_exp_diff))))
# # weight histogram by intensity of both contributing peaks
# hist_weight = hist * np.array([spectrum_exp.peaks.intensities[np.where(mz_exp == mz_exp_i)] * spectrum_exp.peaks.intensities[np.where(mz_exp == mz_exp_i + bin_edges[1])] for mz_exp_i in bin_edges[:-1]]).flatten()
# # plot histogram
# import matplotlib.pyplot as plt
# fig, ax = plt.subplots(1,1)
# plt.bar(bin_edges[:-1], hist_weight, width = 0.5, color='grey', alpha=0.5)
# plt.ylabel('Number of peaks')
# plt.xlabel('m/z difference')
# plt.show()






#spectrum_exp.plot()
# print all scores larger than 0
print([score for score in scores if score["score"] > 0])

scores_times_ln_no_peaks = [score["score"] * np.log(score["matches"]) for score in scores]
indices_rev = np.argsort(scores)[-21:]
#indices_rev = np.argsort(scores_times_ln_no_peaks)[-14:]
indices = indices_rev[::-1]

# print out formula names and shifts of the 5 highest matching spectra
for index in indices:
    print(
        spectra_rel_shifted[index].get("formula"),
        spectra_rel_shifted[index].get("shift"),
        # print out score
        scores[index], # print out score
    )

# get the 5 highest matching spectra
spectra_rel_shifted_top5 = [spectra_rel_shifted[index] for index in indices]

# plot the 5 highest matching spectra







from matplotlib import pyplot as plt
import matplotlib

# intialize plot
fig, ax = plt.subplots(figsize=(1, 1))

#def plot
# plot experimental spectrum as stem plot in black
ax.stem(
    spectrum_exp.peaks.mz,
    spectrum_exp.peaks.intensities,
    label="experimental spectrum",
    linefmt="k-",
    markerfmt="ko",
    basefmt="",
)
# plot spectrum_exp.peaks.mz and -spectrum_exp.peaks.intensities as dark grey crosses
ax.plot(
    spectrum_exp.peaks.mz,
    -spectrum_exp.peaks.intensities,
    "x",
    color="darkgrey",
    label="experimental spectrum",
)
# get first 5 python colors
# colors_user = matplotlib.colors.TABLEAU_COLORS
# plot the 5 highest matching spectra normalized to the experimental spectrum
for i_spec, spectrum_rel_shifted_top5_part in enumerate(spectra_rel_shifted_top5):
    spectrum_rel_shifted_top5_part = normalize_intensities(
        spectrum_rel_shifted_top5_part
    )
    # factor to experiment spectrum:
    # 1) get part of experimental spectrum that contains the masses where the relative spectrum has intensities
    # 2) get the maximum intensity of the part of the experimental spectrum
    # 3) divide the maximum intensity of the part of the experimental spectrum by the maximum intensity of the relative spectrum
    # 4) multiply the relative spectrum by the factor
    # 1)
    # get mz values of experimental spectrum
    mz_exp = spectrum_exp.peaks.mz
    # get mz values of relative spectrum
    mz_rel = spectrum_rel_shifted_top5_part.peaks.mz
    
    spectrum_exp_part = get_part_of_exp_spectrum(mz_exp, mz_rel)

    # get indices of mz values of relative spectrum in mz values of experimental spectrum
    indices_overlapping = np.where(np.isin(mz_exp, mz_rel))
    # if there are no mz values of the relative spectrum in the experimental spectrum, set factor to 0
    if len(indices_overlapping[0]) == 0:
        factor = 0
    elif len(indices_overlapping[0]) == 1:
        # skip this case
        continue
    # if there are mz values of the relative spectrum in the experimental spectrum, set factor to the maximum intensity of the part of the experimental spectrum divided by the maximum intensity of the relative spectrum
    else:
        # get maximum intensity of the part of the experimental spectrum
        max_intensity_exp = spectrum_exp_part.peaks.intensities
        # get maximum intensity of the relative spectrum
        max_intensity_rel = spectrum_rel_shifted_top5_part.peaks.intensities
        # set factor to the maximum intensity of the part of the experimental spectrum divided by the maximum intensity of the relative spectrum
        factor = np.array(max_intensity_exp) / np.array(max_intensity_rel)
        # get weighted average, weight is max_intensity_rel
        factor = np.average(factor, weights=max_intensity_rel)
    # 4)
    # plot the relative spectrum as stem plot
    # get formula name
    formula_name = spectrum_rel_shifted_top5_part.get("formula")
    # mass of the formula
    mass = chem_formula(formula_name).isotope.mass
    # get mass difference to experimental spectrum
    mass_diff =  round(spectrum_rel_shifted_top5_part.peaks.mz[0] - round(mass))

    # get matching score
    score =scores[indices[i_spec]]
    label_name = formula_name + " sc: " + str( round(float(score["score"]), ndigits=3)) + " m: " + str(round(mass)) +  " dm: " + str(mass_diff)
    ax.stem(
        spectrum_rel_shifted_top5_part.peaks.mz,
        np.array(spectrum_rel_shifted_top5_part.peaks.intensities).astype(float)*-factor,
        label=label_name,
        linefmt=['b','g','r','y','m','c','k','b','g','r','y','m','c','k','b','g','r','y','m','c','k','b','g','r','y','m','c','k'][i_spec] + "-",
        markerfmt=['bo','go','ro','yo','mo','co','ko','bv','gv','rv','yv','mv','cv','kv','bp','gp','rp','yp','mp','cp','kp','bD','gD','rD','yD','mD','cD','kD'][i_spec],
        basefmt=" ",
    )
    # add vertical line at mz = 

# add histogram of mz_exp_diff
# ax2 = ax.twinx()
# ax2.bar(bin_edges[:-1], hist, width = 0.5, color='grey', alpha=0.5)
# ax2.set_ylabel('Number of peaks')
# ax2.set_ylim(0, 10)
# # add m/z values below histogram peaks
# for i, v in enumerate(hist):
#     ax2.text(bin_edges[i] - 0.5, v + 0.1, str(int(round(bin_edges[i], ndigits=0))), color='black')





# set x and y limits
ax.set_xlim(0, 300)
ax.set_ylim(-1.05, 1.05)
# set x and y labels
ax.set_xlabel("m/z")
ax.set_ylabel("relative intensity")
# set legend by formula names
# legend_names = [spectrum_rel_shifted_top5_part.get("formula") for spectrum_rel_shifted_top5_part in spectra_rel_shifted_top5]
ax.legend()
# show plot
plt.show()
