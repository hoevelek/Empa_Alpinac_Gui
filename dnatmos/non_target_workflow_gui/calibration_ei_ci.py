"""Script for performing the mass calibration required for the peak extraction scripts.

This one is specially configured for the EI CI files from the 2023 campain.

"""
import logging
from datetime import datetime
import pandas as pd
from alpinac.mass_calibration.mode_mass_calibration import make_mass_calibration, Plots
from alpinac.io_tools_hdf5 import h5_files_from_to
from data_analysis_utils import get_chemical_possible_combinations, get_chemical_possible_combinations_for_elem_list, load_data, get_sum, get_ranges_for_summing, get_calibrant_peaks, mass_index_to_mass, rt_to_rt_index, rt_index_to_rt, mass_to_mass_index
import numpy as np
from matplotlib import pyplot as plt
from molmass import Formula as chem_formula
from scipy.signal import find_peaks, find_peaks_cwt
from pathlib import Path
from GenerateNIST_high_res_file import get_NIST_EI_mass_spectra_by_identifier, load_pseudo_high_res_spectrum
import pubchempy as pcp 
import cirpy
from matchms.importing import load_from_mgf
from alpinac_gui.matchms_funcs import StickSpectra
import shutil

logging.basicConfig()
logging.getLogger('alpinac').setLevel(logging.INFO)


def get_calibrant_formula_masses(calibrant_name:str, calibrant_formula:str, calibrant_rt:float, fragment_masses:float, radicals_allowed:bool=True, path_combinatorial_data:Path=None, RT_calib_params:list = [0.2, 0], dict_targets = {}):
    """
    Get the fragments for the calibrant and append them to the dict_targets dictionary. If the fragments are not yet saved in a file, generate them and save them in a file.
    If radicals are allowed, the fragments are generated with radicals. If not, the fragments are generated as ground-state molecules.
    :param calibrant_name: name of the calibrant
    :param calibrant_formula: formula of the calibrant
    :param calibrant_rt: retention time of the calibrant
    :param fragment_masses: masses of the fragments
    :param radicals_allowed: if True, radicals are allowed in the fragments
    :param path_to_file: path to the file where the precalculated fragments are saved
    :param RT_calib_params: parameters for the RT calibration function: [seconds per bin of RT time, offset of RT time]
    :param dict_targets: dictionary containing the targets
    :return: fragment_formulas: list of the formulas of the fragments, 
             fragment_masses: list of the masses of the fragments, 
             fragments: list of the fragments as strings, 
             calibrant_rt: retention time of the calibrant, 
             calibrant_rt_index: index of the retention time of the calibrant,
             dict_targets: dictionary containing the targets substances
    """
    
    # get element from calibrant formula
    #calibrant_formula = [chem_formula(calibrant_formula).formula]
    elements = [str(atom_comp) for atom_comp in chem_formula(calibrant_formula).composition()]        
    #if path_to_file == None:
    #    path_to_file = Path(path_combinatorial_data) / Path(str("formulas_"+"".join(elements)+"_radicals_"+str(radicals_allowed)+".txt"))
    #try:
    #    data = pd.read_csv(path_to_file, header=None, sep="\t")
    #    fragment_formulas, fragment_masses = data[0].tolist(), data[1].tolist()
    #except:
    _, _, fragment_formulas, fragment_masses = get_chemical_possible_combinations_for_elem_list(elements = elements, target_masses=fragment_masses, radicals_allowed=radicals_allowed, path_to_files=path_combinatorial_data) # TODO: remove restriction to 2 fragments!!!
    #calibrant_formula.extend(fragment_formulas)
    #fragments = [str(formula) for formula in calibrant_formula]
    calibrant_rt_index = rt_to_rt_index(calibrant_rt, RT_calib_params)
    dict_targets[calibrant_name] = [fragment_formulas, calibrant_rt]
    normalized_name=chem_formula(calibrant_formula).formula
    return fragment_formulas, fragment_masses,  normalized_name, calibrant_rt, calibrant_rt_index, dict_targets 

def get_and_plot_target_substances(data, rt_target, rt_cal_params: list = [0.2, 0], mass_cal_params: list = [2835.23579263574, -4049.42781522328], tol=0.1, substance_name = "", fragments = []):
    range_starts, range_ends = get_ranges_for_summing(data, sum_type="RT", threshold_median_multiples=10, thresh_min_range_len= 3)
    print(np.array(range_ends)-np.array(range_starts))
    # eliminiate ranges that are too small
    # print mass spectra for each range
    range_center = 0.5*(range_starts + range_ends)
    line_type_array = ["-", "--", "-.", ":"]*100
    y_max = 0.0
    max_positions = []
    max_values = []

    relevant_masses = []
    subs = fragments
    for sub in subs:
        if isinstance(sub, list):
            #convert to string
            sub = "".join(sub)
        relevant_masses.extend([mass for mass in chem_formula(sub).spectrum()])
    for i, range_start in enumerate(range_starts):
        range_end = range_ends[i]
        rt_sum = get_sum(data, sum_type="RT", range_starts=[range_start], range_ends=[range_end], tol=tol)
        #plt.plot(rt_sum, label="RT range: "+str(range_center[i]))
        #get indeci of rt_sum list
        rt_ind = [i for i in range(len(rt_sum))]
        rt_sum_x, _ = rt_index_to_rt(rt_ind, rt_cal_params)
        # get max y value in relevant range
        rel_indeces = np.where((rt_sum_x > rt_target-10) & (rt_sum_x < rt_target+10))
        y_max_i = np.max(rt_sum[rel_indeces])
        if y_max_i > y_max:
            y_max = y_max_i
        x_max_i = rt_sum_x[rel_indeces][np.argmax(rt_sum[rel_indeces])] 
        # get position of max y value in relevant range
        if(round(mass_index_to_mass(range_center[i], mass_cal_params)[0]) in relevant_masses):
            max_positions.append(x_max_i)
            max_values.append(y_max_i)
            # find max position in interval rt_target-10, rt_target+10 and append to list of max positions using peak_finder
            plt.plot(rt_sum_x, rt_sum, label=round(mass_index_to_mass(range_center[i], mass_cal_params)[0]), linestyle = line_type_array[i])
            # plot vertical line at max position in same color
            plt.axvline(x=x_max_i, linestyle = line_type_array[i], color = "grey")
    # plot horizontal labels
    # if color palette is repeating change the line type
    #plt.legend(bbox_to_anchor=(0.05, 1), loc='upper left', borderaxespad=0.0, fontsize=6, ncol=3)
    # set x axis to rt_ind +- 10
    plt.xlim(rt_target-10, rt_target+10)
    plt.xlabel("RT [s]")
    plt.ylabel("Intensity")
    plt.title("RT spectrum for target substance " + substance_name)
    plt.ylim(0, y_max*1.1)
    plt.legend(ncol=3)
    try:
        subs = fragments
        plt.axvline(x=rt_target, color = "black", linestyle = "--", label = "Estimated RT")
    finally:
        plt.show()
    return max_positions, max_values


#%%
#list all files with ending ".h5" from folder "C:\Users\kaho\Desktop\data\data_Empa\Campaign202303"
path_combinatorial_data = r"C:\Users\kaho\Desktop\data\data_Empa\combinatorial_data"
import os
directory = r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303" # directory to search for files
# mass cal settings
rt_EI = {}
rt_CI = {} 

# for EI look for likely substances, if the are 'not found' (3-parameter fit fails), they will be ignored in calibration 
calibrant_names = ["PFTBA", "Toluene", "Benzene", "Trichlorofluoromethane", "Dichlorodifluoromethane", "Tetrachloromethane"]
#inchi = []
cas_no = {}

cas_no["PFTBA"] = "311-89-7"
cas_no["Toluene"] = "108-88-3"
cas_no["Benzene"] = "71-43-2"
cas_no["Trichlorofluoromethane"] = "75-69-4"
cas_no["Dichlorodifluoromethane"] = "75-71-8"
cas_no["Tetrachloromethane"] = "56-23-5"

#rt_EI['CH2Cl2'] = 1642.1
rt_EI['Benzene'] = 1892.4
rt_EI['Toluene'] = 2213.0
rt_EI["Trichlorofluoromethane"]= 1596.1
rt_EI["Dichlorodifluoromethane"] = 1455.7
rt_EI["Tetrachloromethane"] = 1739.1


directory = r"C:\Users\kaho\Desktop\data\data_Empa\test_mass_cal"
# different mass cal settings 
rt_EI['Benzene'] = 1908
rt_EI['Toluene'] = 2216.0
rt_EI["Trichlorofluoromethane"]= 1597
rt_EI["Dichlorodifluoromethane"] = 1468
rt_EI["Tetrachloromethane"] = 1739.1


directory = r"G:\503_Themen\Klima\ecToFKampagne\Gemessene_Proben\processing_2nd_sample"
# company nearby	date.time.type
# 2	Chemours (Dordrecht)_2	230309.2224.tank
# 4	Givaudon (Kemptal)_2	230310.0855.tank
# 7	Schweizerhalle (Basel)_2	230310.1325.tank
# 8	BASF (Ludwigshafen)_2	230311.0125.tank
# 11	Solvay (Bad Wimpfen)_2	230311.0555.tank
# 14	Eawag ARA (Dübendorf)_2	230315.0350.tank

#1 
rt_EI['Benzene'] = 1908
rt_EI['Toluene'] = 2216.0
rt_EI["Trichlorofluoromethane"]= 1597
rt_EI["Dichlorodifluoromethane"] = 1468
rt_EI["Tetrachloromethane"] = 1739.1

dict_rt_EI_for_file = {}
dict_rt_EI_for_file['230309.2224.tank'] = rt_EI

#4
rt_EI['Benzene'] = 1908.143286
rt_EI['Toluene'] = 2215.6
rt_EI["Trichlorofluoromethane"]= 1597.15
rt_EI["Dichlorodifluoromethane"] = 1467.646519
rt_EI["Tetrachloromethane"] = 1740.09606

dict_rt_EI_for_file['230310.1325.tank'] = rt_EI

#6
rt_EI['Benzene'] = 1908.479381
rt_EI['Toluene'] = 2215.94723
rt_EI["Trichlorofluoromethane"]= 1597.372393
rt_EI["Dichlorodifluoromethane"] = 1468.148516
rt_EI["Tetrachloromethane"] = 1740.769826

dict_rt_EI_for_file['230310.1325.tank'] = rt_EI

#8
rt_EI['Benzene'] = 1908.479381
rt_EI['Toluene'] = 2215.94723
rt_EI["Trichlorofluoromethane"]= 1597.372393
rt_EI["Dichlorodifluoromethane"] = 1468.148516
rt_EI["Tetrachloromethane"] = 1740.769826

dict_rt_EI_for_file['230311.0125.tank'] = rt_EI

#11
rt_EI['Benzene'] = 1908.479381
rt_EI['Toluene'] = 2215.94723
rt_EI["Trichlorofluoromethane"]= 1597.372393
rt_EI["Dichlorodifluoromethane"] = 1468.148516
rt_EI["Tetrachloromethane"] = 1740.769826

dict_rt_EI_for_file['230311.0555.tank'] = rt_EI

#14
rt_EI['Benzene'] = 1908.479381
rt_EI['Toluene'] = 2215.94723
rt_EI["Trichlorofluoromethane"]= 1597.372393
rt_EI["Dichlorodifluoromethane"] = 1468.148516
rt_EI["Tetrachloromethane"] = 1740.769826

dict_rt_EI_for_file['230315.0350.tank'] = rt_EI

# 
dict_rt_EI_for_file['230225.1804.std'] = rt_EI




# first samples

# 1	Chemours (Dordrecht)_1	230308.1924.tank.1 
# 3	Givaudon (Kemptal)_1	230308.2354.tank.3 
# 6	Schweizerhalle (Basel)_1	230309.0554.tank.5 
# 9	BASF (Ludwigshafen)_1	230309.1024.tank.9 
# 12	Solvay (Bad Wimpfen)_1	230309.1624.tank.11 
# 16	Marbach	230310.0254.tank.3


def automatically_extract_calibrant_RT_from_excel(folder, tank):

    subfolder = tank
    filename = tank + '_peak_identification_results_extended.xlsx'
    folder_extended_excel = Path(folder) / Path(subfolder) / Path(filename)

    #load excel to dataframe
    df_calib = pd.read_excel(folder_extended_excel)
    calibrant_names = ['dichloro(difluoro)methane', 'trichloro(fluoro)methane', 'tetrachloromethane', 'benzene', 'toluene']
    # get for each row where of column '0_name' is in calibrant_names the value of column 'RT' and save it in a dictionary
    rt_EI = {}
    for index, row in df_calib.iterrows():
        if row['0_name'] in calibrant_names:
            rt_EI[row['0_name']] = row['RT']    
    # print the dictionary
    for key in rt_EI:
        print(key, rt_EI[key])
    old_names = ['dichloro(difluoro)methane', 'trichloro(fluoro)methane', 'tetrachloromethane', 'benzene', 'toluene']
    new_names = ['Dichlorodifluoromethane', 'Trichlorofluoromethane', 'Tetrachloromethane', 'Benzene', 'Toluene']
    for key_old, key_new in zip(old_names, new_names):
        if key_old in rt_EI:
            rt_EI[key_new] = rt_EI[key_old]
            del rt_EI[key_old]  
    # historical reasons for renaming
    #rt_EI["Trichlorofluoromethane"]= rt_EI['trichloro(fluoro)methane']
    #rt_EI["Dichlorodifluoromethane"] = rt_EI['dichloro(difluoro)methane']
    #rt_EI["Tetrachloromethane"] = rt_EI['tetrachloromethane']
    #rt_EI['Benzene'] = rt_EI['benzene']
    #rt_EI['Toluene'] = rt_EI['toluene']

    # remove old names
    #del rt_EI['trichloro(fluoro)methane']
    #del rt_EI['dichloro(difluoro)methane']
    #del rt_EI['tetrachloromethane']
    #del rt_EI['benzene']
    #del rt_EI['toluene']

    return rt_EI

# 1st sample
directory_sample_files = r"G:\503_Themen\Klima\kaho\preprocessed_data_campaign2023\first_sample_per_location"

# 2nd sample
directory_sample_files = r"G:\503_Themen\Klima\kaho\preprocessed_data_campaign2023\third_sample_per_location"

# 4th sample
directory_sample_files = r"G:\503_Themen\Klima\kaho\preprocessed_data_campaign2023\fourth_sample_per_location"

# 2nd sample
directory_sample_files = r"G:\503_Themen\Klima\kaho\preprocessed_data_campaign2023\second_sample_per_location\Chemours_real_duplicate_2"
file_list_restrictor = "tank"



directory_sample_files = r"G:\503_Themen\Klima\kaho\results_campaign2023\standard2"
auto_extract_RT_from_excel = True
file_list_restrictor = "std"
dict_rt_EI_for_file['230225.1804.std'] = rt_EI
dict_rt_EI_for_file['230403.0617.std'] = rt_EI
dict_rt_EI_for_file['230228.1213.std'] = rt_EI
dict_rt_EI_for_file['230226.1204.std'] = rt_EI
dict_rt_EI_for_file['230225.2234.std'] = rt_EI
rt_EI["tetrachloromethane"] = rt_EI['Tetrachloromethane']




spectra_mode_list = ["EI", "CI"] # CI was done before, but EI failed

if auto_extract_RT_from_excel == True:
    folder_quick_extract = directory_sample_files

    #spectra_mode_list = ["EI", "CI"]
    # create folder old_masscalfiles
    folder_old_masscalfiles = Path(directory_sample_files) / Path("old_masscalfiles")
    if not os.path.exists(folder_old_masscalfiles):
        os.makedirs(folder_old_masscalfiles)
    # move all files ending with _mc.txt to folder old_masscalfiles
    for file in os.listdir(directory_sample_files):
        if file.endswith('EI_mc.txt'):
            try:
                shutil.move(directory_sample_files + "\\" + file, folder_old_masscalfiles)
            except:
                print("File already moved, skipping, delete old_files folder content")


    # list all folders in directory
    folders = [f for f in os.listdir(directory_sample_files) if os.path.isdir(os.path.join(directory_sample_files, f))]
    # remove everythin which does not start with 6 numbers + "."
    tanks = [f for f in folders if f[:6].isdigit()]
    # automatically generate
    for tank in tanks:
        print(tank)
        dict_rt_EI_for_file[tank[:-2]] = automatically_extract_calibrant_RT_from_excel(folder_quick_extract, tank)


# manually overwrite if necessary
# nothing to do here

directory = r"C:\Users\lab134\Desktop\ecToFKampagne\Gemessene_Proben\processing_2nd_sample\queue\Chemours_real_duplicate_2"

file_list = []
ind_for_characteristic_data = 2 # (campain 0,1 was not working (not full RT range))
# loop over each file in the directory
for filename in os.listdir(directory_sample_files):
    if filename.endswith('.h5'):  # check if the file ends with ".py"
        file_list.append(directory_sample_files+"\\"+filename) # print the file name if it meets the condition
        rt_EI = dict_rt_EI_for_file[filename[:-5]]

print(file_list)
# only keep the files that have "tank" in the name
file_list = [file for file in file_list if file_list_restrictor in file]

# Not working: 0, 1, 47

# get calibrant peaks
#%%
DBNIST_WEBPATH = "https://webbook.nist.gov/cgi/cbook.cgi?JCAMP="
download_folder_NIST_pseudo_high_res = r"C:\Users\kaho\Desktop\data\data_Empa\NIST_data\pseudo_high_res_NIST"
download_folder_NIST = r"C:\Users\kaho\Desktop\data\data_Empa\NIST_data\unit_mass_NIST"
path_to_possible_fragments = r"C:\Users\kaho\Desktop\data\data_Empa\combinatorial_data"

# on Bela server

ionisation_dict = {
    # Last value is not included
    "EI": (3,4),
    "CI": (1,2),
}


s_per_bin = float(1/5) # 5Hz
offset = 0 
# typical values for the mass calibration with 2 parameters for approx. plotting
p1 = 2835.23579263574
p2 = -4049.42781522328

# If you want to check for possible calibranr fragments, set plot_summary to True, set ranges if only a part of the file should be used
plot_summary = False
#spectra_mode = "CI"
if plot_summary == True:
    spectra_mode = "CI"
    data, reader = load_data(file_list[ind_for_characteristic_data], ionisation_dict, ionisation_type = spectra_mode)
    mass_sum = get_sum(data, sum_type="mass", tol=0.1)
    fig, ax = plt.subplots(1,1) 
    plt.plot(mass_sum)

    #%%
    RT_sum = get_sum(data, sum_type="RT", tol=0.1)
    # plot the sum
    fig, ax = plt.subplots(1,1)
    plt.plot(RT_sum)


dict_spikes_of_spectra_type = {}
dict_calibrants_EI = {}
dict_calibrants_CI = {}


for identifier in calibrant_names:
    # get inchi keys
    # if not already in the dict:
    if identifier not in cas_no:
        cmps = pcp.get_compounds(identifier, 'name')
        inchis = [cmp.inchi for cmp in cmps]
        cas_no[identifier] = [cirpy.resolve(inchi, 'cas')[-1] for inchi in inchis] # cactus server (http://cactus.nci.nih.gov/chemical/structure/) might be down, was down since more than a week


# test if the rt times are to far away, by getting them from the excel sheet
#


# sort rt_EI by rt
calibrant_sorted_by_RT = sorted(rt_EI, key=rt_EI.get)

# PLOT AND CHECK THE CALIBRANT PEAKS: to figure out best rts
#fig, ax = plt.subplots(1,1)
#max_pos_Toluol, max_val_Toluol = get_and_plot_target_substances(data, rt_Toluol, s_per_bin, offset, p1, p2, tol=0.1, substance_name="Toluol", fragments  = fragments_Toluol)

# GET PSEUDO-HIGH_Res SPECTRA FROM NIST PLUS CONVERSION
# get cas
# https://webbook.nist.gov/cgi/cbook.cgi?JCAMP=C311897&Index=0&Type=Mass 

# make sure to copy the files in the path_hr folder. For heavy molecules it might takes some minutes to calculate
for key in cas_no:
    #try to find high_res file
    cas = cas_no[key]
    path_save = str(Path(download_folder_NIST_pseudo_high_res) / Path(cas + ".mgf"))
    pseudo_high_res, pseudo_high_res_largest_N, formulas_largest_N = load_pseudo_high_res_spectrum(cas, path_save, largest_N_peaks=12)
    fig, ax = plt.subplots()
    pseudo_high_res_largest_N.plot_stem(ax, label_peaks = "formulas")
    #set title to iupac name and cas no
    plt.title(key + " " + cas)
    directory_mass_cal = r"C:\Users\kaho\Desktop\data\data_Empa\mass_cal\calibrant_spectra\EI"
    # save figure
    plt.savefig(os.path.join(directory_mass_cal, f"{key}_{cas_no[key]}.png"))

    dict_calibrants_EI[key] = formulas_largest_N
    # write into file named key_+cas_no+.txt into folder  the formula fragments
    with open(os.path.join(directory_mass_cal, f"{key}_{cas_no[key]}.txt"), 'w') as f:
        for formula in formulas_largest_N:
            f.write(f"{formula}\n")
    

# In case of PFTBA some of the peaks are ?wrongly? assigned to sth containing N atoms (not distinguishable isotopic patterns), overwrite:
# see https://link.springer.com/chapter/10.1007/978-3-319-54398-7_3
dict_calibrants_EI["PFTBA"] = [
 'CF3',
# 'C3FN',
 'C2F4',
# 'C4F2N',
# 'C2F4N',
 'C2F5',
 'C3F5',
# 'C5F3N',
 'C4F9',
# 'C6F7N',
# 'C5F10N',
 'C8F16N', #added, as at a good mass (but small)
# 'C9F20N']
]

dict_calibrants_CI["PFTBA"] = [
 'CF3',
 'C2F4',
 #'C4F2N',
 'C2F5',
 'C3F5',
 'C4F9',
 'C8F16N', #large in CI
]

# SET THE TARGETS FOR THE CALIBRATION FOR EI

# SET THE TARGETS FOR THE CALIBRATION FOR CI
dict_targets_CI = {}
dict_targets_CI['0.0'] = dict_calibrants_CI["PFTBA"] # continous
dict_spikes_of_spectra_type['CI'] = dict_targets_CI


file_list = ['C:\\Users\\kaho\\Desktop\\data\\data_Empa\\Campaign202303\\230224.1630.air.1.h5']
# remove object from filelist for which already a calibration file exists
file_list = [file for file in file_list if not os.path.isfile(file[:-3]+"mc.txt")]

for file in file_list:
    for spectra_mode in spectra_mode_list:
        print("processing", file)
        print("spectra_mode", spectra_mode)

        if spectra_mode == "EI":
            dict_spikes = {}

            # GET CALIBRANT RETENTION TIMES BY FINDING PEAKS IN THE RT SPECTRUM AT THE START AND END OF THE RT SPECTRUM
            # The spikes do not come at the same times, so one has to search them for each file
            data, reader = load_data(file, ionisation_dict, ionisation_type = spectra_mode)
            # plot RT spectrum
            #RT_sum = get_sum(data, sum_type="RT", tol=0.1)
            # plot the sum
            #fig, ax = plt.subplots(1,1)
            #plt.plot(RT_sum)

            _, _, _, calibrant_peaks_Rts_start, spikes_xFWHM_left, spikes_xFWHM_right = get_calibrant_peaks(data, start_or_end="start", factor_median_for_threshold=10, plot=False)
            # get the time from index
            spike_durations_start = abs(np.array(spikes_xFWHM_left) - np.array(spikes_xFWHM_right))
            _, _, _, calibrant_peaks_Rts_end, spikes_xFWHM_left, spikes_xFWHM_right = get_calibrant_peaks(data, start_or_end="end", factor_median_for_threshold=10, plot=False)
            #FWHM as spike duration, much smaller than the real RT peak duration (can not use for calibration function RT tolerance)
            spike_durations_end = abs(np.array(spikes_xFWHM_left) - np.array(spikes_xFWHM_right))

            # #add the calibrant peaks to the dict_spikes
            for i in range(len(calibrant_peaks_Rts_start)):
                dict_spikes[str(round(calibrant_peaks_Rts_start[i], ndigits=2))] = dict_calibrants_EI["PFTBA"]

            for cal in calibrant_sorted_by_RT:
                str_rt = str(round(rt_EI[cal], ndigits=2)) # to be processed by mygu script
                dict_spikes[str_rt] = dict_calibrants_EI[cal]

            for i in range(len(calibrant_peaks_Rts_end)):
                dict_spikes[str(round(calibrant_peaks_Rts_end[i], ndigits=2))] = dict_calibrants_EI["PFTBA"]
            
            dict_spikes_of_spectra_type['EI']  = dict_spikes
            tof_corr = None
        else:
            tof_corr = lambda x: 1.001 * x

        try:
            make_mass_calibration(
                path_file=Path(file),
                #segments=(1,2),
                segments=ionisation_dict[spectra_mode],
                mass_cal_mode_index=2,
                dict_spikes=dict_spikes_of_spectra_type[spectra_mode],
                spike_duration=10,
                mass_diff_ppm_threshold = 40,
                #averaging_timspan=600, # Longer timspan
                #spacing_timespan=400,
                mass_domain_range=0.4,
                tofid_base_correction=tof_corr,
                plots=[Plots.MASS_DRIFT, Plots.PARAMETERS],
                show_plots=False,
                spectra_mode=spectra_mode
            ) 
        except:
            print("Error in file", file)
    # generate the excel file

