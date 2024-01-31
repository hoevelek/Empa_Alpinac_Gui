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
import os

logging.basicConfig()
logging.getLogger('alpinac').setLevel(logging.INFO)

# Add calibrant peaks to the mass calibration folder

# Main Code

dict_rt_EI_for_file = {}
rt_EI = {}
cas_no = {}
# calibrant_names = ["PFTBA", "Toluene", "Benzene", "Trichlorofluoromethane", "Dichlorodifluoromethane", "Tetrachloromethane"]
cas_no["PFTBA"] = "311-89-7"
cas_no["Toluene"] = "108-88-3"
cas_no["Benzene"] = "71-43-2"
cas_no["Trichlorofluoromethane"] = "75-69-4"
cas_no["Dichlorodifluoromethane"] = "75-71-8"
cas_no["Tetrachloromethane"] = "56-23-5"
cas_no["PFTBA"] = "311-89-7"
calibrant_names = cas_no.keys()
dict_calibrants_EI = {}
dict_calibrants_CI = {}

# Organize calibrant RT values
#folder_quick_extract = directory_sample_files
#folder_extended_excel = Path(folder_quick_extract) / Path("tank" + '_peak_identification_results_extended.xlsx')
#dict_rt_EI_for_file['tank'] = automatically_extract_calibrant_RT_from_excel(folder_quick_extract, "tank")

# Organize calibrant peaks
for identifier in calibrant_names:
    if identifier not in cas_no:
        cmps = pcp.get_compounds(identifier, 'name')
        inchis = [cmp.inchi for cmp in cmps]
        cas_no[identifier] = [cirpy.resolve(inchi, 'cas')[-1] for inchi in inchis]

# Load and plot pseudo-high-res spectra
download_folder_NIST_pseudo_high_res = r"C:\Users\kaho\Desktop\data\data_Empa\NIST_data\pseudo_high_res_NIST"
for key in cas_no:
    cas = cas_no[key]
    path_save = str(Path(download_folder_NIST_pseudo_high_res) / Path(cas + ".mgf"))
    pseudo_high_res, pseudo_high_res_largest_N, formulas_largest_N = load_pseudo_high_res_spectrum(cas, path_save, largest_N_peaks=12)
    fig, ax = plt.subplots()
    pseudo_high_res_largest_N.plot_stem(ax, label_peaks="formulas")
    plt.title(key + " " + cas)
    directory_mass_cal = r"C:\Users\kaho\Desktop\data\data_Empa\mass_cal\calibrant_spectra\EI"
    plt.savefig(os.path.join(directory_mass_cal, f"{key}_{cas_no[key]}.png"))
    dict_calibrants_EI[key] = formulas_largest_N
    with open(os.path.join(directory_mass_cal, f"{key}_{cas_no[key]}.txt"), 'w') as f:
        for formula in formulas_largest_N:
            f.write(f"{formula}\n")

# MANUALLY ADD CALIBRANTS IF DIFFERENT FROM LARGEST 12 NIST PEAKS
# write file same file for calibrants with are special and do not correspond to normal NIST spectra, e.g. CI etc
# key = "PFTBA"
# dict_calibrants_EI["PFTBA"] = [
#  'CF3',
# # 'C3FN',
#  'C2F4',
# # 'C4F2N',
# # 'C2F4N',
#  'C2F5',
#  'C3F5',
# # 'C5F3N',
#  'C4F9',
# # 'C6F7N',
# # 'C5F10N',
#  'C8F16N', #added, as at a good mass (but small)
# # 'C9F20N']
# ]
# # save to file in same scheme and folder as EI
# with open(os.path.join(directory_mass_cal, f"{key}_{cas_no[key]}.txt"), 'w') as f:
#     for formula in dict_calibrants_EI["PFTBA"]:
#         f.write(f"{formula}\n")

directory_mass_cal = r"C:\Users\kaho\Desktop\data\data_Empa\mass_cal\calibrant_spectra\CI"

dict_calibrants_CI["PFTBA"] = [
 'CF3',
 'C2F4',
 #'C4F2N',
 'C2F5',
 'C3F5',
 'C4F9',
 'C8F16N', #large in CI
]
with open(os.path.join(directory_mass_cal, f"{key}_{cas_no[key]}.txt"), 'w') as f:
    for formula in dict_calibrants_CI["PFTBA"]:
        f.write(f"{formula}\n")