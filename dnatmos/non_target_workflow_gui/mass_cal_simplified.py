import logging
from datetime import datetime
import pandas as pd
from alpinac.mass_calibration.mode_mass_calibration import make_mass_calibration, Plots
# import extracting mass peaks routines
from alpinac.mode_extraction 
from alpinac.io_tools_hdf5 import h5_files_from_to
from data_analysis_utils import get_chemical_possible_combinations, get_chemical_possible_combinations_for_elem_list, load_data, get_sum, get_ranges_for_summing, get_calibrant_peaks, mass_index_to_mass, rt_to_rt_index, rt_index_to_rt, mass_to_mass_index, get_RT_spec_for_mass_range, get_mass_spec_for_RT_range
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import axes
from molmass import Formula as chem_formula
from scipy.signal import find_peaks, find_peaks_cwt
from pathlib import Path
from GenerateNIST_high_res_file import get_NIST_EI_mass_spectra_by_identifier, load_pseudo_high_res_spectrum
import pubchempy as pcp 
import cirpy
from matchms.importing import load_from_mgf
from alpinac_gui.matchms_funcs import StickSpectra
import shutil
from pathlib import Path
from alpinac_processing_utils import create_modified_input_file, get_alpinac_results, get_matchms_results_of_alpinac_fragments

# set up logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Define a dictionary to store loaded data
loaded_data_cache = {}


DBNIST_WEBPATH = "https://webbook.nist.gov/cgi/cbook.cgi?JCAMP="



# user paths and defaults
user_paths = {}
user_paths["raw_data_file"] = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230224.1630.air.1.h5")
user_paths["mass_cal"] = Path(r"C:\Users\kaho\Desktop\data\data_Empa\mass_cal\")
user_paths["calibrant_spectra"] = Path(r"C:\Users\kaho\Desktop\data\data_Empa\mass_cal\calibrant_spectra")
user_paths["databases_empa"] = Path(r"C:\Users\kaho\Desktop\data\data_Empa\databases_Empa\TargetList") # target list & spectra must be in folder!!
user_paths["data_bases_extra_CH"] = Path(r"C:\Users\kaho\Desktop\data\data_Empa\databases_Empa\TargetList\nist_autodownloads\renamed_files")
user_paths["alpinac_input_files"] = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230224.1630.air.1\alpinac_input_files")
user_paths["alpinac_output_files"] = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230224.1630.air.1\alpinac_output_files")
user_paths["combinatorial_data"] = Path(r"C:\Users\kaho\Desktop\data\data_Empa\combinatorial_data")
user_paths["download_NIST"] = Path(r"C:\Users\kaho\Desktop\data\data_Empa\NIST_data\unit_mass_NIST")
user_paths["download_NIST_pseudo_high_res"] = Path(r"C:\Users\kaho\Desktop\data\data_Empa\NIST_data\pseudo_high_res_NIST")


defaults = {}
defaults_mass_cal = {}
defaults_mass_cal["mass_ranges_start"] = []
defaults_mass_cal["mass_ranges_end"] = []
defaults_mass_cal["rt_index_to_rt_mode"] = [0.]
defaults_mass_cal["rt_index_to_rt_params"] = [float(1/5),0]
defaults_mass_cal["mass_index_to_mass_mode"] = [2]
defaults_mass_cal["mass_index_to_mass_params"] = [2835.23579263574, -4049.42781522328, 0.5]
defaults_mass_cal["range_tol"] = [0.1]
defaults_mass_cal["normalize"] = [0]
defaults_mass_cal["mass_for_rt_dependent_calibration"] = [75]  # as calibration is RT, mass-dependent for the plotting only a mean is taken
defaults_mass_cal["rt_for_mass_dependent_calibration"] = [1000]  # as calibration is RT, mass-dependent for the plotting only a mean is taken

defaults_extraction = {}
defaults_extraction["rt_start_extract"] = [1000]
defaults_extraction["rt_stop_extract"] = [2000]
defaults_extraction["adduct"] = [None]
defaults_extraction["mass_list"] = [None]
defaults_extraction["ionization_segments"] = [(1,2)]
defaults_extraction["rt_sigma_guess_scale"] = [1.0]
defaults_extraction["mass_u_ppm_min"] = [0]

defaults_alpinac = {}
defaults_alpinac["target_elements"] = ["C", "H", "O", "N", "S", "Cl", "Br", "F", "P"]
defaults_alpinac["fitting_method"] = ["continuous"]


def get_content_of_calibrant_spectra(file_path):
    # read the file
    file_df = pd.read_csv(file_path, names=["fragments"])
    #convert to list
    file_list = file_df['fragments'].tolist()
    return file_list


def file_content_rt_of_suited_substances(directory_path, file_name: "000_default.txt"):
    # if files exists, load it
    if os.path.isfile(directory_path + file_name):
        file_df = pd.read_csv(directory_path + file_name, sep=" = ", names=["substance", "RT"])
    else:
        file_df = pd.read_csv(directory_path + "000_default.txt", sep=" = ", names=["substance", "RT"])
        #warining: file not found
        logger.warning("Set to default as file not found: " + directory_path + file_name)
    return file_df


def load_data_bases(db_dir, db_dir_extra_CH):
    # Load the data bases
    dbs: dict[str, dict[str, Spectrum]] = {
        "nist": JdxReader(db_dir / "nist_db").read_all(),
        "myriam": AlpinacFragmentsReader(db_dir / "myriam_ei").read_all(),
        "xtra_CH": JdxReader(db_dir_extra_CH).read_all(),
    }
    # target_list (RT)
    use_substances_with_unknown_RT = False
    df_peaks = pd.read_csv(db_dir / "TargetList.csv")
    df_pseudo_peaks = pd.read_csv(db_dir / "PseudoTargetList.csv")
    if use_substances_with_unknown_RT:
        df_peaks = pd.concat([df_peaks, df_pseudo_peaks], ignore_index=True)


def load_data_ext(file_path, ionisation_dict, ionisation_type="EI"):
    """
    Load data from a file.

    Parameters:
        file_path (str): Path to the data file.
        ionisation_dict (dict): Dictionary containing ionisation ranges.
        ionisation_type (str): Type of ionisation (default is "EI").

    Returns:
        tuple: A tuple containing the loaded data and the reader object.
    """
    # Check if data is already in the cache
    if file_path in loaded_data_cache:
        return loaded_data_cache[file_path]

    # If not in the cache, load the data
    data, reader = load_data_internal(file_path, ionisation_dict, ionisation_type)

    # Store the loaded data in the cache
    loaded_data_cache[file_path] = (data, reader)

    return data, reader


def load_data_internal(file_path, ionisation_dict, ionisation_type):
    # Replace this with your actual data loading logic
    data, reader = load_data(raw_data_file, ionisation_dict, ionisation_type = ionization_mode)

    return data, reader


def show_data_raw_RT_plot(raw_data_file, defaults, ionisation_dict, ionisation_type="EI", mass_ranges_start: list = [], mass_ranges_end: list = [], ax = None):
    """
    Show the raw data.

    Parameters:
        raw_data_file (str): Path to the data file.
        defaults (dict): Dictionary containing the default parameters.
        ionisation_dict (dict): Dictionary containing ionisation ranges.
        ionisation_type (str): Type of ionisation (default is "EI").
        mass_ranges_start (list): List of mass ranges to sum over.
        mass_ranges_end (list): List of mass ranges to sum over.
        ax (matplotlib.axes.Axes): Axes object to plot into (default is None).
    """
    flag = False
    # load data
    data, reader = load_data_ext(raw_data_file, ionisation_dict, ionisation_type = ionization_mode)

    # defaults
    rt_index_to_rt_mode = int(defaults["rt_index_to_rt_mode"][0])
    rt_index_to_rt_params = defaults["rt_index_to_rt_params"]
    mass_index_to_mass_mode = defaults["mass_index_to_mass_mode"][0]
    mass_index_to_mass_params = defaults["mass_index_to_mass_params"]
    range_tol = defaults["range_tol"][0]
    normalize = bool(defaults["normalize"][0])
    if mass_ranges_start == []:
        mass_ranges_start = None
    if mass_ranges_end == []:
        mass_ranges_end = None

    #if in reader, use reader values
    try:
        mean_mass = np.mean(np.array(mass_ranges_start) + 0.5*(np.array(mass_ranges_end) - np.array(mass_ranges_start)))
    except:
        mean_mass = defaults["mass_for_rt_dependent_calibration"][0]

    try:
        mass_index_to_mass_mode = reader.mass_calibration_data[ionisation_type].mass_cal_mode
        mass_index_to_mass_params_help = reader.mass_calibration_data[ionisation_type].get_parameters_for_rt_nearest_neighbours(70) 
        mass_index_to_mass_params = [mass_index_to_mass_params_help[key] for key in mass_index_to_mass_params_help.keys() if key.startswith("p")]
    except:
        pass


    # if mass_ranges_start is [] set to None
    if mass_ranges_start == []:
        mass_ranges_start = None
    if mass_ranges_end == []:
        mass_ranges_end = None
    # use get_sum to get the sum of all scans
    result = get_RT_spec_for_mass_range(
        data=data,
        ranges_start=mass_ranges_start,
        ranges_end=mass_ranges_end,
        rt_index_to_rt_mode=rt_index_to_rt_mode,
        rt_index_to_rt_params=rt_index_to_rt_params,
        mass_index_to_mass_mode=int(mass_index_to_mass_mode),
        mass_index_to_mass_params=mass_index_to_mass_params,
        range_tol=range_tol,
        range_mode='abs',
        sum_x='RT',
        normalize=False,
        )
    # Get the time-of-flight axis
    try:
        RTime = reader.get_rt_axis()
    except:
        RTime = result["RT"]  # works of calibration did not work
    # Get the mass axis
    # plot the sum
    if ax is None:
        fig, ax = plt.subplots(1,1)
        flag = True
    if normalize:
        ax.plot(RTime, result["intensity"]/np.max(result["intensity"]))
    else:
        ax.plot(RTime, result["intensity"])
    ax.set_xlabel("Retention time (s)")
    ax.set_ylabel("Intensity")
    ax.set_title("Raw data")
    if flag:
        plt.show()

    return RTime, result["intensity"]


def show_data_raw_mass_plot(raw_data_file, defaults, ionisation_dict, ionisation_type="EI", RT_ranges_start: list = [], RT_ranges_end: list = [], ax = None):
    """
    Show the raw data along the mass axis.

    Parameters:
        raw_data_file (str): Path to the data file.
        defaults (dict): Dictionary containing the default parameters.
        ionisation_dict (dict): Dictionary containing ionisation ranges.
        ionisation_type (str): Type of ionisation (default is "EI").
        mass_ranges_start (list): List of mass ranges to sum over.
        mass_ranges_end (list): List of mass ranges to sum over.
        ax (matplotlib.axes.Axes): Axes object to plot into (default is None).
    """
    # load data
    data, reader = load_data_ext(raw_data_file, ionisation_dict, ionisation_type=ionization_mode)

    # defaults
    rt_index_to_rt_mode = int(defaults["rt_index_to_rt_mode"][0])
    rt_index_to_rt_params = defaults["rt_index_to_rt_params"]
    mass_index_to_mass_mode = int(defaults["mass_index_to_mass_mode"][0])
    mass_index_to_mass_params = defaults["mass_index_to_mass_params"]
    range_tol = defaults["range_tol"][0]
    normalize = bool(defaults["normalize"][0])
    if RT_ranges_start == []:
        RT_ranges_start = None
    if RT_ranges_end == []:
        RT_ranges_end = None

    # calculate the mean RT
    try:
        mean_RT = np.mean(np.array(RT_ranges_start) + 0.5*(np.array(RT_ranges_end) - np.array(RT_ranges_start)))
    except:
        mean_RT = defaults["rt_for_mass_dependent_calibration"][0]
    
    # mass axis
    try:
        mass = reader.get_mass_axis_of(mean_RT)
    except:
        # use calibration
        mass = mass_index_to_mass(
            mass_index=np.arange(data.shape[1]),
            mode=mass_index_to_mass_mode,
            params=mass_index_to_mass_params,
        )

    # get the sum
    result = get_mass_spec_for_RT_range(
        data=data,
        ranges_start=RT_ranges_start,
        ranges_end=RT_ranges_end,
        rt_index_to_rt_mode=rt_index_to_rt_mode,
        rt_index_to_rt_params=rt_index_to_rt_params,
        mass_index_to_mass_mode=mass_index_to_mass_mode,
        mass_index_to_mass_params=mass_index_to_mass_params,
        range_tol=range_tol,
        range_mode='abs',
        sum_x='mass',
        normalize=False,
        )

    # show the plot
    if ax is None:
        fig, ax = plt.subplots(1,1)
        flag = True
    if normalize:
        ax.plot(mass, result["intensity"]/np.max(result["intensity"]))
    else:
        ax.plot(mass, result["intensity"])
    ax.set_xlabel("Mass (u)")
    ax.set_ylabel("Intensity")
    ax.set_title("Raw data, ")
    if flag:
        plt.show()

    return mass, result["intensity"]


def generate_calibration_dict(path, calibrant_injected, calib_mode, ionization_mode, mass_range):
    """
    Generates the dict_spikes and dict_spikes_mass for the mass calibration
    """
    dict_spikes = {}
    dict_spikes_of_spectra_type = {}


    for file in os.listdir(path):
        if file.startswith(calibrant_injected) and file.endswith(".txt"):
            calibrant_injected_file = file

    if calib_mode == "continous":
        # read from file
        dict_spikes['0.0'] = get_content_of_calibrant_spectra(path + calibrant_injected_file)
        tof_corr = lambda x: 1.001 * x
        dict_spikes_of_spectra_type[ionization_mode]  = dict_spikes

    elif calib_mode == "spiked":
            # GET CALIBRANT RETENTION TIMES BY FINDING PEAKS IN THE RT SPECTRUM AT THE START AND END OF THE RT SPECTRUM
            # The spikes do not come at the same times, so one has to search them for each file
            data, reader = load_data_ext(raw_data_file, ionisation_dict, ionisation_type = ionization_mode)
            # plot RT spectrum
            #RT_sum = get_sum(data, sum_type="RT", tol=0.1)
            # plot the sum
            #fig, ax = plt.subplots(1,1)
            #plt.plot(RT_sum)
            calibrant_substances_df = file_content_rt_of_suited_substances(directory_path, Path(raw_data_file).name.replace(".h5", "calibrant."+ionization_mode+".txt"))
            calibrant_sorted_by_RT = calibrant_substances_df.sort_values(by=['RT'])

            _, _, _, calibrant_peaks_Rts_start, spikes_xFWHM_left, spikes_xFWHM_right = get_calibrant_peaks(data, start_or_end="start", factor_median_for_threshold=10, plot=True)
            # get the time from index
            spike_durations_start = abs(np.array(spikes_xFWHM_left) - np.array(spikes_xFWHM_right))
            _, _, _, calibrant_peaks_Rts_end, spikes_xFWHM_left, spikes_xFWHM_right = get_calibrant_peaks(data, start_or_end="end", factor_median_for_threshold=10, plot=False)
            #FWHM as spike duration, much smaller than the real RT peak duration (can not use for calibration function RT tolerance)
            spike_durations_end = abs(np.array(spikes_xFWHM_left) - np.array(spikes_xFWHM_right))

            # #add the calibrant peaks to the dict_spikes
            for i in range(len(calibrant_peaks_Rts_start)):
                dict_spikes[str(round(calibrant_peaks_Rts_start[i], ndigits=2))] = get_content_of_calibrant_spectra(path + calibrant_injected_file)

            for ind, (cal, RT) in calibrant_sorted_by_RT.iterrows():
                print(cal)
                print(RT)
                str_rt = str(round(RT, ndigits=2))
                if str_rt not in dict_spikes.keys():
                    cal_file = ""
                    for file in os.listdir(path):
                        if file.startswith(cal) and file.endswith(".txt"):
                            cal_file = file
                    dict_spikes[str_rt] = get_content_of_calibrant_spectra(path + cal_file)

            for i in range(len(calibrant_peaks_Rts_end)):
                dict_spikes[str(round(calibrant_peaks_Rts_end[i], ndigits=2))] = get_content_of_calibrant_spectra(path + calibrant_injected_file)
            tof_corr = None
            dict_spikes_of_spectra_type[ionization_mode]  = dict_spikes
    else: 
        logger.error("calib_mode not known, use 'continous' or 'spiked'")
    return dict_spikes_of_spectra_type


    # calculate the mass of all fragments in the dict_spikes and remove those that are not in the mass range
    dict_spikes_mass = {}
    for key in dict_spikes.keys():
        dict_spikes_mass[key] = []
        for fragment in dict_spikes[key]:
            mass = chem_formula(fragment).isotope.mass
            if mass >= mass_range[0] and mass <= mass_range[1]:
                dict_spikes_mass[key].append(mass)
            else:
                print("Fragment not in mass range: ", fragment, mass)
                # remove the fragment from the list
                dict_spikes[key].remove(fragment)
    

def do_mass_calibration(raw_data_file, ionization_mode, calib_mode, calibrant_injected, mass_range):

    
    """
    Does the mass calibration for one file
    """
    # generate the dict_spikes and dict_spikes_mass
    dict_spikes_of_spectra_type = generate_calibration_dict(path, calibrant_injected, calib_mode, ionization_mode, mass_range)
    # do the mass calibration

    try:
        make_mass_calibration(
            path_file=Path(raw_data_file),
            #segments=(1,2),
            segments=ionisation_dict[ionization_mode],
            mass_cal_mode_index=2,
            dict_spikes=dict_spikes_of_spectra_type[ionization_mode],
            spike_duration=10,
            mass_diff_ppm_threshold = 40,
            #averaging_timspan=600, # Longer timspan
            #spacing_timespan=400,
            mass_domain_range=0.4,
            tofid_base_correction=tof_corr,
            plots=[Plots.MASS_DRIFT, Plots.PARAMETERS],
            show_plots=False,
            spectra_mode=ionization_mode,
        ) 
    except:
        print("Error in file", raw_data_file)


def do_extraction(data_file_directory, rt_start_extract, rt_stop_extract, outfile_path, outfile_name, adduct = None, mass_list=None, ionization_segments=(1,2), rt_sigma_guess_scale=1.0, mass_u_ppm_min=0):
    alp_input_file = input_file_folder / Path(alp_input_filename)
    try:
        if alp_input_file.exists():
            print("input file already exists, skipping")
        else:
            print("alpinac input file does not exist, creating")
            # h5path = r'c:\Users\kaho\Desktop\data\data_Empa\Campaign202303_interesting_file\230311.0125.tank.1.h5'
            # elution_time = 1655
            # input_file_folder = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303_interesting_file\230311.0125.tank.1\alpinac_input_files")
            # alp_input_filename = 'test_cmp12.txt'
            # adduct = 'H'
            # rt_range_left = 2.0
            # rt_range_right = 3.0
            # mass_list = None
            nontarget_peak_extraction(
                path_file=h5path,
                rt_start_extract=elution_time - rt_range_left,
                rt_stop_extract=elution_time + rt_range_right,
                adduct=adduct,
                outfile_path=input_file_folder,
                outfile_name=alp_input_filename,
                mass_list=mass_list,
                ionization_segments=ionization_segments,
                rt_sigma_guess_scale=1.0,  # might give wrong result for large RT
                mass_u_ppm_min=0,
            )
            print("input file {} generated".format(alp_input_file))
    except:
        print("Error in extraction process")
        # set to failed


def do_alpinac(alp_input_file, path_for_alpinac_output, target_elements, fitting_method):
    """
    Runs alpinac on one file
    """
    try:
        path_for_alpinac_output_cmp = Path(path_for_alpinac_output)/Path((str(alpinac_input_file.name)[:-4] + "_mod").replace("input", "output"))
        if not path_for_alpinac_output_cmp.exists():
            path_for_alpinac_output_cmp.mkdir()
                                                                                
        make_identification_from_file(path_file=alpinac_input_file_mod, output_dir=path_for_alpinac_output_cmp, show_plt_figures=False, target_elements=target_elements, fitting_method=FittingMethod.continuous, timeout_time=30*60)
        print("alpinac ran successfully")
        # set alpinac_status to 1 for cmp_id
    except:
        print("alpinac failed")
    

if __name__ == "__main__":

    # GUI is split into 3 panels: upper, middle, lower
    # upper panel is split into 4 panels: left (called raw data), middle left (called calibration), middle right (called extraction), right (called identification)

    # Upper task bar: "Paths", "Calibrate Defaults", "Extract Defaults", "Identify Defaults": 
    # if click on "Paths": open a window with name fields for all paths and a dropdown for each path to select the path
    # if click on "Calibrate Defaults": open a window with name fields for all calibrate defaults and let the user enter the values, save on click on "save and close" and close window
    # if click on "Extract Defaults": open a window with name fields for all extract defaults and let the user enter the values, save on click on "save and close" and close window
    # if click on "Identify Defaults": open a window with name fields for all identify defaults and let the user enter the values, save on click on "save and close" and close window
    # will save to file

    # defaults will be used if not overwritten by user input in the GUI
    # Step 0.0: Load defaults from file (if file does not exist, use the defaults defined in the code)

    # Upper panel
    # Step 1: Select raw data file 
    # Step 2: Fill "EI segments" and "CI segments" field with tuple: e.g (1,2) or (1)
    # Step 3: Call show_data_raw_RT_plot for EI and CI in a common plot (y1 and y2 axis) and show it in the full middle panel. (One should be able to zoom in a plot)
        # right top of the plot two fields: "mass start" and "mass stop" (default: blank, blank)
        # as soon as the user enters a value, the plot is updated

        # E.g. if both fileds blank for EI
        show_data_raw_RT_plot(raw_data_file, defaults, ionisation_dict, ionization_mode, mass_ranges_start=None, mass_ranges_end=None, ax=None)
        # E.g. if values in both fields for EI
        show_data_raw_RT_plot(raw_data_file, defaults, ionisation_dict, ionization_mode, mass_ranges_start=np.array([30]), mass_ranges_end=np.array([70]), ax=None)

    # Step 4: Call show_data_raw_mass_plot for EI and CI in a common plot (y1 and y2 axis) and show it in the full lower panel at the same time as the RT plot.
        # right top of the plot two fields: "RT start" and "RT stop" (default: blank, blank)
        # as soon as the user enters a value, the plot is updated (One should be able to zoom in a plot)
        # E.g. if both fileds blank for EI
        show_data_raw_mass_plot(raw_data_file, defaults, ionisation_dict, ionization_mode, RT_ranges_start=None, RT_ranges_end=None, ax=None)
        # E.g. if values in both fields for EI
        show_data_raw_mass_plot(raw_data_file, defaults, ionisation_dict, ionization_mode, RT_ranges_start=np.array([1000]), RT_ranges_end=np.array([2000]), ax=None)

    # Step 5: If user presses "Calibrate" button, call do_mass_calibration
    # Next to calibrate button: "Calibration mode" (continous, spiked)
    # Next to calibrate button: "Calibrant injected" (PFTBA) (drop down)
    # Next to calibrate button: Use known substances (dropdown to folder with known substances)
    # Next to calibrate button: "Message field" (e.g. "Calibration successful")

    # Step 6: If user presses "Extract" button, call do_extraction
    # Next to extract button: "RT start" (default: blank)
    # Next to extract button: "RT stop" (default: blank)
    # Next to extract button: "Adduct" (default: blank)
    # Next to extract button: "Mass list" (default: blank)
    # Next to extract button: "Message field" (e.g. "Extraction successful")

    # Step 7: If user presses "Identify Alpinac" button, call do_alpinac with the input file created in step 6
    # if sucessful call get_alpinac_results
    # Next to identify button: "Target elements" (default: C, H, O, N, S, Cl, Br, F, P)
    # Next to identify button: "Fitting method" (default: continuous)
    # Next to identify button: "Message field" (e.g. "Alpinac failed" or if sucessful return of get_alpinac_results)

    # Step 8: If user presses "Identify Matching" button, call do_matchms with the input file created in step 6
    # if sucessful call get_matchms_results_of_alpinac_fragments
    # Next to identify button: "Message field" (e.g. "Matching failed" or if sucessful return of get_matchms_results_of_alpinac_fragments)


    





#directory_path = "C:\\Users\\kaho\\Desktop\\data\\data_Empa\\mass_cal\\"


#folder_calibrants = "C:\\Users\\kaho\\Desktop\\data\\data_Empa\\mass_cal\\calibrant_spectra\\"
#folder_rt_calibrants = "C:\\Users\\kaho\\Desktop\\data\\data_Empa\\mass_cal\\"
#ionisation_dict = {"EI": (3,4), "CI": (1,2)}    
#ionisation_dict = {"EI": (3,4)}
# EI only:
#ionization_mode = "EI"
#calib_mode = "spiked"
#calibrant_injected = "PFTBA" # either spiked or continous


file_list = ['C:\\Users\\kaho\\Desktop\\data\\data_Empa\\Campaign202303\\230224.1630.air.1.h5']
file_list = ['C:\\Users\\kaho\\Desktop\\data\\data_Empa\\Campaign202303\\230331.0646.tank.1.h5']

raw_data_file = file_list[0]

# 4 cases: EI continous, EI spiked, CI continous, CI spiked
# In case of EI continous, the calibrant is injected continously, so dict_targets['0.0'] = calibrant
path = directory_path + "calibrant_spectra\\" + ionization_mode + "\\"
mass_range = (20, 300)


#download_folder_NIST_pseudo_high_res = r"C:\Users\kaho\Desktop\data\data_Empa\NIST_data\pseudo_high_res_NIST"
#download_folder_NIST = r"C:\Users\kaho\Desktop\data\data_Empa\NIST_data\unit_mass_NIST"
#path_to_possible_fragments = r"C:\Users\kaho\Desktop\data\data_Empa\combinatorial_data"

#db_dir = Path(r"C:\Users\kaho\Desktop\data\data_Empa\databases_Empa\TargetList")
#db_dir_extra_CH = Path(r"C:\Users\kaho\Desktop\data\data_Empa\databases_Empa\TargetList\nist_autodownloads\renamed_files")