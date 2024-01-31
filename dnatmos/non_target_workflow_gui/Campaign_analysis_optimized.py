## This file should do optimized analysis of the data from the campaign. It will work on single lines of the extracted data.

## % How to proceed
# 1. do mass calibration for EI and CI using calibration_ei_ci.py
# 2. detect peaks, extract the and compare with Ei data with library entries, write results into excel file
# 3. use interesting compounds, (no carbohydrates, less then 90% confidelity interval) and run alpinac on them
# 4. write results into excel file

#%% Load excel file with library entries
import ast
from functools import partial
import logging
from multiprocessing import managers
import threading
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
import shutil
import multiprocessing


# mass calib with external progr.
# List all folders in directory "C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\" and read the files in them

import os
import glob
from alpinac.mass_calibration.utils import IonizationMethod
from alpinac.mode_extraction.mode_extraction_nontarget import nontarget_peak_extraction
from alpinac.mode_identification.main_from_file import make_identification_from_file
from alpinac.mode_identification.compound_identification import compound_identification
from alpinac.utils_data_extraction import FittingMethod
import alpinac_sups.read_h5 as rh5
from alpinac_gui.matchms_funcs import formula_from_namestring, formula_to_latex, halocarbname2formula, plot_mass_spectra, load_from_mgf, metadata_processing, peak_processing, AlpinacData, halocarbname2formula
import pubchempy as pcp
import cirpy
from molmass import Formula as chem_formula
from pyvalem.formula import Formula
import chardet
from pathlib import Path
import multiprocessing as mp # paralellized processing
from multiprocessing import Manager, Pool


# create extended table
from Campaign_quick_extract import prepare_extended_data_table, formula_by_name
from data_analysis_utils import formula_by_name, clean_empa_name_string
        
from functools import partial
import multiprocessing

from alpinac_processing_utils import create_modified_input_file, get_alpinac_results, get_matchms_results_of_alpinac_fragments

logging.basicConfig(level=logging.INFO)
logging.getLogger('alpinac').setLevel(logging.INFO)



def do_extraction(extended_quick_res_path, ionization_segments, rt_range_left:float = 2.0, rt_range_right:float = 3.0, adduct: str = 'H', cmp_ids: list = None, mass_list = None):
    # loaf extended excel file
    if not isinstance(cmp_ids, np.ndarray):
        cmp_ids = np.array(cmp_ids)
    # check if file already exists
    if os.path.isfile(extended_quick_res_path):
        df_ext = pd.read_excel(Path(extended_quick_res_path), sheet_name="Sheet1")
    else:
        raise Exception("extended excel file does not exist")
    # get hdf5 file path
    h5path = Path(extended_quick_res_path).parent.parent/ (Path(extended_quick_res_path).parent.name + ".h5")

    # if not already exist create folder for alpinac input files
    if not os.path.exists(h5path):
        # warning
        print("h5 file not found at: " + str(h5path))
    else:
        # if folder alpinac_input_files does not exist, create it
        input_file_folder = Path(h5path).parent/Path(Path(h5path).name).stem / "alpinac_input_files"
        if not os.path.exists(input_file_folder):
            os.mkdir(input_file_folder)
            print("alpinac_input_files folder {} created".format(input_file_folder))

        # do extraction for all cmps for which extraction_status == 0
        print(cmp_ids)
        if cmp_ids[0] == None:
            cmp_ids = df_ext.loc[df_ext["extraction_status"] == 0, "cmp_id"].to_numpy()

        # do extraction for all cmps for which extraction_status == 0
        for cmp_id in cmp_ids:
            # do extraction
            print("doing extraction for cmp_id: " + str(cmp_id))
            # get index of cmp_id in df
            ind_cmp = df_ext.index[df_ext["cmp_id"] == cmp_id].to_numpy()[0]
            # get RT of cmp_id
            elution_time = df_ext.loc[ind_cmp, "RT"]
            # call extraction function
            alp_input_filename = df_ext.loc[ind_cmp, "alpinac_input_file_EI_CI"]
            alp_input_file = input_file_folder/Path(alp_input_filename)
            try:
                df_ext = pd.read_excel(Path(extended_quick_res_path), sheet_name="Sheet1")
                print(df_ext)
                if alp_input_file.exists():
                    print("input file already exists, skipping")
                    #continue
                else:
                    print("alpinac input file does not exist, creating")
                    nontarget_peak_extraction(
                        path_file = h5path,                              
                        rt_start_extract = elution_time-rt_range_left,                              
                        rt_stop_extract = elution_time+rt_range_right, 
                        adduct= adduct, 
                        outfile_path=input_file_folder,
                        outfile_name=alp_input_filename,
                        #plot_pseudovoigt = True,
                        #coelution_plot_mode = True,
                        #mass_bin_plot_mode = True,
                        #save_plots=[Plots.CHROMATOGRAPHY_FITTING, Plots.COELUTION_PLOT, Plots.RT_PEAK_DETECTION],
                        #mass_start = mass_start,                              
                        #mass_stop = mass_stop,                              
                        mass_list = mass_list,
                        ionization_segments=ionization_segments,
                        rt_sigma_guess_scale=1.0, # might give wrong result for large RT 
                        mass_u_ppm_min = 10.0,
                    )
                    print("input file {} generated".format(alp_input_file))
                    df_ext.loc[ind_cmp, "extraction_status"] = 1
            except:
                print("Error in extraction process")
                #set to failed
                df_ext.loc[ind_cmp, "extraction_status"] = 2
            df_ext.to_excel(extended_quick_res_path, index=False)
    return df_ext

def extract_cmp(cmp_id, extended_quick_res_path, ionization_segments, rt_range_left, rt_range_right, adduct, h5path, input_file_folder, df_ext, mass_list):
    # do extraction
    print("doing extraction for cmp_id: " + str(cmp_id))
    # get index of cmp_id in df
    ind_cmp = df_ext.index[df_ext["cmp_id"] == cmp_id].to_numpy()[0]
    # get RT of cmp_id
    elution_time = df_ext.loc[ind_cmp, "RT"]
    # call extraction function
    alp_input_filename = df_ext.loc[ind_cmp, "alpinac_input_file_EI_CI"]
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
        df_ext.loc[ind_cmp, "extraction_status"] = 1
    except:
        print("Error in extraction process")
        # set to failed
        df_ext.loc[ind_cmp, "extraction_status"] = 2

def do_extraction_parallel(extended_quick_res_path, ionization_segments, rt_range_left=2.0, rt_range_right=3.0, adduct='H', cmp_ids=None, mass_list=None):
    # load extended excel file
    if not isinstance(cmp_ids, np.ndarray):
        cmp_ids = np.array(cmp_ids)
    # check if file already exists
    if os.path.isfile(extended_quick_res_path):
        df_ext = pd.read_excel(Path(extended_quick_res_path), sheet_name="Sheet1")
    else:
        raise Exception("extended excel file does not exist")

    # get hdf5 file path
    h5path = Path(extended_quick_res_path).parent.parent / (Path(extended_quick_res_path).parent.name + ".h5")

    # if not already exist create folder for alpinac input files
    if not os.path.exists(h5path):
        # warning
        print("h5 file not found at: " + str(h5path))
    else:
        # if folder alpinac_input_files does not exist, create it
        input_file_folder = Path(h5path).parent / Path(Path(h5path).name).stem / "alpinac_input_files"
        if not os.path.exists(input_file_folder):
            os.mkdir(input_file_folder)
            print("alpinac_input_files folder {} created".format(input_file_folder))

        # do extraction for all cmps for which extraction_status == 0
        print(cmp_ids)
        if cmp_ids[0] is None:
            cmp_ids = df_ext.loc[df_ext["extraction_status"] == 0, "cmp_id"].to_numpy()

        # Create a pool of worker processes
        pool = mp.Pool()

        # Define the function arguments for parallel execution
        args = [(cmp_id, extended_quick_res_path, ionization_segments, rt_range_left, rt_range_right, adduct, h5path, input_file_folder, df_ext, mass_list) for cmp_id in cmp_ids]

        # Map the worker function to the arguments using the pool
        pool.starmap(extract_cmp, args)

        # Close the pool of worker processes and wait for them to complete
        pool.close()
        pool.join()

        # Save the updated dataframe to the extended_quick_res_path
        #df_ext.to_excel(extended_quick_res_path, index=False)

    return df_ext




def run_alpinac(alpinac_input_file, path_for_alpinac_output_list, target_element_list, alpinac_mode_list):
    # read the alpinac input file name
    alpinac_input_file_mod = Path(str(alpinac_input_file)[:-4] + "_mod.txt")

    # run alpinac for each mode
    #for path_for_alpinac_output, target_elements, alpinac_mode in zip(path_for_alpinac_output_list, target_element_list, alpinac_mode_list):
    for path_for_alpinac_output, target_elements in zip(path_for_alpinac_output_list, target_element_list):


        try:
            path_for_alpinac_output_cmp = Path(path_for_alpinac_output)/Path((str(alpinac_input_file.name)[:-4] + "_mod").replace("input", "output"))
            if not path_for_alpinac_output_cmp.exists():
                path_for_alpinac_output_cmp.mkdir()
                                                                                 
            make_identification_from_file(path_file=alpinac_input_file_mod, output_dir=path_for_alpinac_output_cmp, show_plt_figures=False, target_elements=target_elements, fitting_method=FittingMethod.continuous, timeout_time=30*60)
            print("alpinac ran successfully")
            # set alpinac_status to 1 for cmp_id
        except:
            print("alpinac failed")
            # set alpinac_status to 2 for cmp_id

def run_alpinac_slow(cmp_id, file, df_alp, path_for_alpinac_output_list, target_element_list, alpinac_mode_list):
    # read the alpinac input file name
    alpinac_input_file = Path(file).parent / "alpinac_input_files" / df_alp.loc[df_alp["cmp_id"] == cmp_id, "alpinac_input_file_EI_CI"].values[0]
    alpinac_input_file_mod = Path(str(alpinac_input_file)[:-4] + "_mod.txt")

    # run alpinac for each mode
    for path_for_alpinac_output, target_elements, alpinac_mode in zip(path_for_alpinac_output_list, target_element_list, alpinac_mode_list):
        try:
            path_for_alpinac_output_cmp = Path(path_for_alpinac_output)/Path((str(alpinac_input_file.name)[:-4] + "_mod").replace("input", "output"))
            if not path_for_alpinac_output_cmp.exists():
                path_for_alpinac_output_cmp.mkdir()
                                                                                 
            make_identification_from_file(path_file=alpinac_input_file_mod, output_dir=path_for_alpinac_output_cmp, show_plt_figures=False, target_elements=target_elements)
            print("alpinac ran successfully")
            # set alpinac_status to 1 for cmp_id
            #df_alp.loc[df_alp["cmp_id"] == cmp_id, alpinac_mode] = 1
            #alp_res_str = get_alpinac_results(path_for_alpinac_output)
            #df_alp.loc[df_alp["cmp_id"] == cmp_id, alpinac_mode] = alp_res_str
        except:
            print("alpinac failed")
            # set alpinac_status to 2 for cmp_id
            #df_alp.loc[df_alp["cmp_id"] == cmp_id, alpinac_mode] = 2

def parallel_alpinac(file):
    df_ext = pd.read_excel(Path(file), sheet_name="Sheet1", index_col=0)
    

    # check for which cmp_ids the extraction_status is 1 and modified input file exists
    cmp_ids_extracted_mod = df_ext.loc[(df_ext["extraction_status"] == 1), "cmp_id"].to_numpy()    
    cmp_ids_not_yet_eval = df_ext.loc[df_ext["alpinac_status"] == 0, "cmp_id"].to_numpy()
    cmp_ids = np.intersect1d(cmp_ids_extracted_mod, cmp_ids_not_yet_eval)

    # generate list of aplianac input files as dict
    alpinac_input_files_dict = {}
    for cmp_id in cmp_ids:
        alpinac_input_files_dict[cmp_id] = Path(file).parent / "alpinac_input_files" / df_ext.loc[df_alp["cmp_id"] == cmp_id, "alpinac_input_file_EI_CI"].values[0]

    
 

    # initialize alpinac modes
    path_for_alpinac_output_list = [Path(file).parent / "alpinac_results", Path(file).parent / "alpinac_results_no_N"]
    target_element_list = [None, "CHOBrClFIPS"]
    alpinac_mode_list = ["alpinac_EI_CI", "alpinac_EI_CI_no_N"]

    # if not already exist create folders for alpinac OUTPUT files, one for EI_CI and one for EI_CI_no_N
    for path in path_for_alpinac_output_list:
        if not os.path.exists(path):
            os.mkdir(path)

    # Create a pool of worker processes
    pool = mp.Pool()

    # Define the function arguments for parallel execution
    args = [(alpinac_input_files_dict[cmp_id], path_for_alpinac_output_list, target_element_list, alpinac_mode_list) for cmp_id in cmp_ids]

    # Map the worker function to the arguments using the pool
    pool.starmap(run_alpinac, args)

    # Close the pool of worker processes and wait for them to complete
    pool.close()
    pool.join()

    # overwrite excel file with new results
    
    

def timeout_handler(signum, frame):
    # This function will be called when the timeout is reached
    raise TimeoutError("Function timed out")

def process_cmp(alpinac_input_file_mod, cmp_id, path_for_alpinac_output, target_elements):
    # Get input file name
    # Run alpinac
    path_for_alpinac_output_dir = Path(str(path_for_alpinac_output / (Path(alpinac_input_file_mod.name).stem).replace("input", "output")))
    status = 0
    try:
        # Call the function that you want to set a timeout for
        make_identification_from_file(path_file=alpinac_input_file_mod, output_dir=path_for_alpinac_output_dir, show_plt_figures=False, target_elements=target_elements, fitting_method=FittingMethod.discrete, timeout_time=30*30) 
        # Set alpinac_status to 1 for cmp_id
        print("Alpinac ran successfully")    
        status = 1  
    except TimeoutError as e:
        # Handle the timeout exception here
        status = 4
        print(e)
    # except any other exception
    except:
        print("Alpinac failed")
        # Set alpinac_status to 2 for cmp_id
        status = 2
    finally:
        # Cancel the timeout alarm
        pass

    return status

def process_with_timeout(partial_func, timeout, alpinac_input_file_mod, cmp_id, path_for_alpinac_output, target_elements):
    def target():
        nonlocal result
        result = partial_func(alpinac_input_file_mod, cmp_id, path_for_alpinac_output, target_elements)

    result = None
    thread = threading.Thread(target=target)
    thread.start()
    thread.join(timeout=timeout)

    if thread.is_alive():
        thread.join()  # Wait for the thread to complete (optional)
        raise TimeoutError(f"Process exceeded {timeout} seconds.")
    return result




if __name__ == '__main__':
    logging.basicConfig(level=logging.WARNING) 
    logger = logging.getLogger("dnatmos_analysis")

    use_substances_with_unknown_RT = True

    mp.freeze_support()

    # on Laptop
    # directory for database
    #db_dir = Path(r"G:\503_Themen\Klima\TargetList")
    db_dir = Path(r"C:\Users\kaho\Desktop\data\data_Empa\databases_Empa\TargetList")
    #db_dir_extra_CH = Path(r"G:\503_Themen\Klima\TargetList\nist_autodownloads\renamed_files")
    db_dir_extra_CH = Path(r"C:\Users\kaho\Desktop\data\data_Empa\databases_Empa\TargetList\nist_autodownloads\renamed_files")
    folder_path_campaign = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303_second_twin")

    # on calculation server (Bela Rechneré G): "broken pipe, maybe caused by Network problems"
    db_dir = Path(r"G:\503_Themen\Klima\ecToFKampagne\Gemessene_Proben\databases_Empa\TargetList")
    db_dir_extra_CH = Path(r"G:\503_Themen\Klima\ecToFKampagne\Gemessene_Proben\databases_Empa\TargetList\nist_autodownloads\renamed_files")
    folder_path_campaign = Path(r"G:\503_Themen\Klima\ecToFKampagne\Gemessene_Proben\processing_2nd_sample")

    # copied to Desktop
    db_dir = Path(r"c:\Users\lab134\Desktop\ecToFKampagne\Gemessene_Proben\databases_Empa\TargetList")
    db_dir_extra_CH = Path(r"c:\Users\lab134\Desktop\ecToFKampagne\Gemessene_Proben\databases_Empa\TargetList\nist_autodownloads\renamed_files")
    folder_path_campaign = Path(r"c:\Users\lab134\Desktop\ecToFKampagne\Gemessene_Proben\processing_2nd_sample\queue")

    path_target_list = db_dir/Path("TargetList_extended.csv")

    folder_path_campaign = Path(r"c:\Users\lab134\Desktop\ecToFKampagne\Gemessene_Proben\processing_4th_sample\queue")
    folder_path_campaign = Path(r"c:\Users\lab134\Desktop\ecToFKampagne\Gemessene_Proben\processing_1st_sample\queue") # remove input files, mass calibration change
    #folder_path_campaign = Path(r"c:\Users\lab134\Desktop\ecToFKampagne\Gemessene_Proben\processing_3rd_sample\queue")
    folder_path_campaign = Path(r"c:\Users\lab134\Desktop\ecToFKampagne\Gemessene_Proben\processing_2nd_sample\queue\Chemours_real_duplicate_2") # remove input files, mass calibration change

    # on server
    
    # get path for those files having quick extraction passes, if not do quick extraction first, they will not show up here
    #folder_path_campaign = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303_interesting_file")
    folder_path_campaign = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303_second_twin") # discrete fitting
    folder_path_campaign = Path(r"G:\503_Themen\Klima\kaho\results_campaign2023\processed")
    folder_path_campaign = Path(r"G:\503_Themen\Klima\kaho\results_campaign2023\Chemours_real_duplicate_2")
    folder_path_campaign = Path(r"G:\503_Themen\Klima\kaho\results_campaign2023\Schweizerhalle_real_duplicate_2")
    folder_path_campaign = Path(r"G:\503_Themen\Klima\kaho\results_campaign2023\BASF_real_duplicate_2")
    folder_path_campaign = Path(r"G:\503_Themen\Klima\kaho\results_campaign2023\Solvay_real_duplicate_2")

    # machine learning calibration, single
    folder_path_campaign = Path(r"G:\503_Themen\Klima\kaho\results_campaign2023\BASF_single_1")
    folder_path_campaign = Path(r"G:\503_Themen\Klima\kaho\results_campaign2023\Schweizerhalle_single_1")
    folder_path_campaign = Path(r"G:\503_Themen\Klima\kaho\results_campaign2023\Solvay_single_2")
    folder_path_campaign = Path(r"G:\503_Themen\Klima\kaho\results_campaign2023\Givaudon_single_1") # was not completely calculated, redo on server!
    folder_path_campaign = Path(r"G:\503_Themen\Klima\kaho\results_campaign2023\Chemours_real_duplicate_2")

    folder_path_campaign = Path(r"c:\Users\lab134\Desktop\ecToFKampagne\Gemessene_Proben\standard_good_mass_cal")




    # Load the data bases
    dbs: dict[str, dict[str, Spectrum]] = {
        "nist": JdxReader(db_dir / "nist_db").read_all(),
        "myriam": AlpinacFragmentsReader(db_dir / "myriam_ei").read_all(),
        "xtra_CH": JdxReader(db_dir_extra_CH).read_all(),
    }
    # target_list (RT)
    df_peaks = pd.read_csv(db_dir / "TargetList.csv")
    df_pseudo_peaks = pd.read_csv(db_dir / "PseudoTargetList.csv")
    if use_substances_with_unknown_RT:
        df_peaks = pd.concat([df_peaks, df_pseudo_peaks], ignore_index=True)


    df_peaks_extended = pd.read_csv(path_target_list)
    # add Nan columns with iupac_name and inchi_key tp df_pseudo_peaks
    df_pseudo_peaks_extended = df_pseudo_peaks.copy()
    df_pseudo_peaks_extended["iupac_name"] = np.nan
    df_pseudo_peaks_extended["inchi_key"] = np.nan
    df_peaks_extended = pd.concat([df_peaks_extended, df_pseudo_peaks_extended], ignore_index=True)


    

    # dictionary for names of library which have not-standard names:

    # DICTIONARIES # NOT NICE, COMMON ONE FOR QUICK EXTRACT FILE AND THIS ONE
    # dictionary with manually added formulas
    formulas_dict_manually_added = {
        'alpha-methylstyrene-b':['C2Cl3F3', 'prop-1-en-2-ylbenzene', '7407', 'InChI=1S/C9H10/c1-8(2)9-6-4-3-5-7-9/h3-7H,1-2H3', '2095.7'],
        'dichloro-trifluoro-pyridine':['C5Cl2F3N', '3,5-dichloro-2,4,6-trifluoropyridine', '61273', 'InChI=1S/C5Cl2F3N/c6-1-3(8)2(7)5(10)11-4(1)9', '2745.6'],
        'C5H6O-a':['C5H6O', 'None', '10103117','None', '1856.7'],
        'C5H6O-b':['C5H6O', '2-methylfuran', '10797', 'InChI=1S/C5H6O/c1-5-3-2-4-6-5/h2-4H,1H3', '1865.9'],
        'DMS':['C2H6S', 'methylsulfanylmethane', '1068', 'InChI=1S/C2H6S/c1-3-2/h1-2H3', '1572'],
        'ethylbenz':['C8H10', 'ethylbenzene', '7500', 'InChI=1S/C8H10/c1-2-8-6-4-3-5-7-8/h3-7H,2H2,1H3','2776.7'],
        'perfluorodimethyl-cyclohexane':['C8F16', '1,1,2,2,3,3,4,5,5,6-decafluoro-4,6-bis(trifluoromethyl)cyclohexane', '78975', 'InChI=1S/C8F16/c9-1(7(19,20)21)3(11,12)2(10,8(22,23)24)5(15,16)6(17,18)4(1,13)14]', '1779.1'],
        'perfluorodecaline':['C10F18', '1,1,2,2,3,3,4,4,4a,5,5,6,6,7,7,8,8,8a-octadecafluoronaphthalene', '9386', 'InChI=1S/C10F18/c11-1-2(12,5(17,18)9(25,26)7(21,22)3(1,13)14)6(19,20)10(27,28)8(23,24)4(1,15)16]', '2073.2'],
        'i-C6F14pent':['C6F14', 'None', 'None', 'InChI=1S/C6F14/c7-1(4(12,13)14,5(15,16)17)2(8,9)3(10,11)6(18,19)20', '1639.5'],
    }


    # GLOBAL VARIABLES
    formula_dict = {}
    
    try:
        with open(db_dir.parent /Path("target_list_supplementary_files")/ 'formula_dict.txt', 'r') as f:
            for line in f:
                # if line contains one ":"
                if line.count(':') == 1:
                    (key, val) = line.split(':')
                    print(key)
                    print(val)
                    # covert string to list
                    val = val.strip().strip('[]').split('\', ')
                    # convert list elements to string
                    val = [x.strip().strip("'") for x in val]
                    # assign the list elements to the variables
                    formula_dict[key] = val
    except:
        print('No formula dictionary found, creating one')


    # get path for those files having quick extraction passes, if not do quick extraction first, they will not show up here
    #folder_path_campaign = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303_interesting_file")
    folder_path_campaign = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303_second_twin") # discrete fitting
    folder_path_campaign = Path(r"G:\503_Themen\Klima\kaho\results_campaign2023\processed")
    folder_path_campaign = Path(r"G:\503_Themen\Klima\kaho\results_campaign2023\Chemours_real_duplicate_2")
    folder_path_campaign = Path(r"G:\503_Themen\Klima\kaho\results_campaign2023\Schweizerhalle_real_duplicate_2")
    folder_path_campaign = Path(r"G:\503_Themen\Klima\kaho\results_campaign2023\BASF_real_duplicate_2")
    folder_path_campaign = Path(r"G:\503_Themen\Klima\kaho\results_campaign2023\Solvay_real_duplicate_2")

    # machine learning calibration, single
    folder_path_campaign = Path(r"G:\503_Themen\Klima\kaho\results_campaign2023\BASF_single_1")
    folder_path_campaign = Path(r"G:\503_Themen\Klima\kaho\results_campaign2023\Schweizerhalle_single_1")
    folder_path_campaign = Path(r"G:\503_Themen\Klima\kaho\results_campaign2023\Solvay_single_2")
    folder_path_campaign = Path(r"G:\503_Themen\Klima\kaho\results_campaign2023\Givaudon_single_1") # was not completely calculated, redo on server!
    folder_path_campaign = Path(r"G:\503_Themen\Klima\kaho\results_campaign2023\Chemours_real_duplicate_2")

    
    dir_contents = os.listdir(folder_path_campaign)
    # only those starting with a date
    dir_contents = list(filter(lambda x: x.startswith("2"), dir_contents))
    all_excel_files = []
    all_excel_files_extended = []
    extension1 = "*results.xlsx"
    extension2 = "*results_extended.xlsx"
    status_dict = {"alpinac_EI_CI": "alpinac_status", "alpinac_EI_CI_no_N": "alpinac_no_N_status"} # alpinac status

    # define ionization segments
    ionization_segments={
        IonizationMethod.EI: (3,4),
        IonizationMethod.CI: (1,2),
    }

    # iterate over each item in the directory
    for item in dir_contents:
        # construct the full path to the item
        item_path = os.path.join(folder_path_campaign, item)

        # check if the item is a folder
        if os.path.isdir(item_path):
            print("processing folder: " + item)
            print(item_path)
            for extension, list_excel in zip([extension1, extension2], [all_excel_files, all_excel_files_extended]):
                # construct the full search path using the folder path and extension
                search_path = os.path.join(item_path, extension)
                # use the glob module to find all files that match the search path
                excel_files = glob.glob(search_path)
                #check if excelfile starts with a date in format 231231
                if excel_files:
                    excel_files = [excel_file for excel_file in excel_files if Path(excel_file).name[:6].isdigit()]
                # extend if not empty
                if excel_files:
                    list_excel.extend(excel_files)

    do_extract = False
    if do_extract == True:
        # do data extraction for all excel files
        for file_ext in all_excel_files_extended:
            # check for which cmp_ids the extraction_status is 0
            df_ext = pd.read_excel(Path(file_ext), sheet_name="Sheet1")

            # clear mass_calibration data
            for ionization_mode in ['EI', 'CI']:
                mass_calib_file = Path(str(Path(file_ext).parent) + "." +ionization_mode + "_mc" + ".txt")
                if os.path.exists(mass_calib_file):
                    # save copy as backup
                    mass_calib_bu = str(mass_calib_file)[:-4] + "_backup.txt"
                    shutil.copy(mass_calib_file, mass_calib_bu)

                    # delete all conten from original file
                    #with open(mass_calib_file, "w") as f:
                    #    pass
                    
                    counter_rm = 0
                    # remove lines which contain "[]"
                    with open(mass_calib_file, "r") as f:
                        lines = f.readlines()
                    with open(mass_calib_file, "w") as f:
                        for line in lines:
                            if "[]" not in line:
                                f.write(line)
                            else:
                                counter_rm += 1
                    print("removed {} lines from {}".format(counter_rm, mass_calib_file))
                    #open file again and replace first digit with digit - counter_rm
                    with open(mass_calib_file, "r") as f:
                        lines = f.readlines()
                    with open(mass_calib_file, "w") as f:
                        for i, line in enumerate(lines):
                            if i == 0:
                                ind_end_digit = line.find("\t")
                                f.write(str(int(line[:ind_end_digit]) - counter_rm) + line[ind_end_digit:])
                            else:
                                f.write(line)
                    print("removed {} RTs contributing to mass cal of file {}".format(counter_rm, mass_calib_file))

        
            cmp_ids = df_ext.loc[df_ext["extraction_status"] == 0, "cmp_id"].to_numpy()

            # do extraction
            #df_ext = do_extraction(extended_quick_res_path=file_ext, ionization_segments=ionization_segments, rt_range_left = 2.0, rt_range_right= 3.0, adduct= 'H', cmp_ids=cmp_ids, mass_list = None)
            df_ext = do_extraction_parallel(extended_quick_res_path=file_ext, ionization_segments=ionization_segments, rt_range_left = 2.0, rt_range_right= 3.0, adduct= 'H', cmp_ids=cmp_ids, mass_list = None)
    
    update_excel_file = True
    if update_excel_file == True:
        for file_ext in all_excel_files_extended:
            # check for which cmp_ids the extraction_status is 1
            df_ext = pd.read_excel(Path(file_ext))
            alp_input_filenames = [df_ext.loc[df_ext["cmp_id"] == cmp_id, "alpinac_input_file_EI_CI"].values[0] for cmp_id in df_ext["cmp_id"]]
            extraction_sucess = [1 if (Path(file_ext).parent / "alpinac_input_files"/ alp_input_filename).exists() else 0 for alp_input_filename in alp_input_filenames]
            # set and save extraction_status
            df_ext["extraction_status"] = extraction_sucess
            #save df_ext to excel file
            df_ext.to_excel(file_ext, index=False)

    create_mod_input = True
    if create_mod_input == True:
        # create modified input files for all excel files
        for file_ext in all_excel_files_extended:
            # check for which cmp_ids the extraction_status is 1
            df_ext = pd.read_excel(Path(file_ext), sheet_name="Sheet1")
            # for each cmp check if alp_input_file_EI_CI exists and if yes set extraction_status to 1
            cmp_ids = df_ext.loc[df_ext["extraction_status"] == 1, "cmp_id"].to_numpy()
            # create modified input files
            for cmp_id in cmp_ids:
                # read alpinac input file name
                alp_input_file_relative = df_ext.loc[df_ext["cmp_id"] == cmp_id, "alpinac_input_file_EI_CI"].values[0]
                # create modified input file
                df_ext = create_modified_input_file(
                    alp_input_filename=Path(file_ext).parent / "alpinac_input_files"/ alp_input_file_relative,
                    extended_quick_res_path=Path(file_ext),
                    df_ext=df_ext,
                    cmp_id=cmp_id,
                )
                df_ext.to_excel(file_ext, index=False)


    # run_alpinac_parallel_without_timeout = False
    # if run_alpinac_parallel_without_timeout == True:

    #     # Run alpinac for all excel files
    #     for file in all_excel_files_extended:
    #         # Read the alpinac input file name
    #         df_ext = pd.read_excel(Path(file), sheet_name="Sheet1", index_col=0)
    #         df_alp = df_ext.copy()
    #         status_list = df_alp["alpinac_status"].to_list()
            
    #         cmp_ids_not_yet_processed = df_alp.loc[df_alp["alpinac_status"] == 0, "cmp_id"].to_numpy()
    #         cmp_ids_extracted_mod = df_alp.loc[(df_alp["extraction_status"] == 1), "cmp_id"].to_numpy()
    #         cmp_ids = np.intersect1d(cmp_ids_extracted_mod, cmp_ids_not_yet_processed)

    #         # Initialize alpinac modes
    #         alpinac_mode_list = ["alpinac_EI_CI", "alpinac_EI_CI_no_N"]

    #         # dictionaries for alpinac results
    #         path_for_alpinac_output_dict = {"alpinac_EI_CI": Path(file).parent / "alpinac_results", "alpinac_EI_CI_no_N": Path(file).parent / "alpinac_results_no_N"}
    #         target_elements_dict = {"alpinac_EI_CI": None, "alpinac_EI_CI_no_N": "CHOBrClFIPS"}
            
    #         # Create folders for alpinac OUTPUT files if they don't exist
    #         for path in path_for_alpinac_output_dict:
    #             if not path_for_alpinac_output_dict[path].exists():
    #                 path_for_alpinac_output_dict[path].mkdir()

    #         # generate alpinac_input_file_mod_dict:
    #         alpinac_input_file_mod_dict = {}
    #         for cmp_id in cmp_ids:
    #             # read alpinac input file name
    #             alp_input_file_relative = df_alp.loc[df_alp["cmp_id"] == cmp_id, "alpinac_input_file_EI_CI"].values[0]
    #             # create modified input file
    #             alpinac_input_file_mod_dict[cmp_id]=Path(file).parent / "alpinac_input_files"/ Path(alp_input_file_relative.replace(".txt", "_mod.txt"))

    #         # Create a pool of worker processes
    #         # loop over all 
    #         for alpinac_mode in alpinac_mode_list:
    #             # Create a pool of worker processes
    #             with Pool() as pool:
    #                 # Create a shared df
    #                 manager = Manager()
    #                 shared_df = manager.list(status_list) # Create a shared DataFrame using Manager
    #                 # Define the partial function with shared_df as an argument
    #                 partial_process_cmp = partial(process_cmp, shared_df=shared_df)

    #                 # Map the partial function to the cmp_ids using multiple processes
    #                 pool.starmap(partial_process_cmp, [(alpinac_input_file_mod_dict[cmp_id], cmp_id, path_for_alpinac_output_dict[alpinac_mode], target_elements_dict[alpinac_mode]) for cmp_id in cmp_ids])
    #                 pool.close()
    #                 pool.join()

    #             # Save the final results
    #             df_alp[status_dict[alpinac_mode]] = shared_df
    #             df_alp.to_excel(file[:-5] + "_extended.xlsx", index=False)

        


    update_alpinac_status = True
    if update_alpinac_status == True:
        alpinac_mode = ["alpinac_EI_CI", "alpinac_EI_CI_no_N"]
        # load excel files
        for file in all_excel_files_extended:
            path_for_alpinac_output_dict = { "alpinac_EI_CI": Path(file).parent / "alpinac_results", "alpinac_EI_CI_no_N": Path(file).parent / "alpinac_results_no_N"}
            df_ext = pd.read_excel(Path(file))
            #get status 
            for alp_md in alpinac_mode:
                # Read the alpinac input file name
                # check if at least one folder named starting with "Compound" exists in resultsfolder
                # Create folders for alpinac OUTPUT files if they don't exist
                if path_for_alpinac_output_dict[alp_md].exists():
                    # go through all folders in path_for_alpinac_output_dict[path]
                    for folder in path_for_alpinac_output_dict[alp_md].iterdir():
                        # get cmp_id
                        print(folder)
                        cmp_id = folder.name.replace("_mod","").strip().split("_")[-1]
                        # list subfolders
                        subfolders = [subfolder for subfolder in folder.iterdir() if subfolder.is_dir()]
                        # check if at least one folder named starting with "Compound" exists in resultsfolder
                        if any([subfolder.name.startswith("Compound") for subfolder in subfolders]):
                            #check if subfolder is not empty
                            if any([len(os.listdir(subfolder)) > 0 for subfolder in subfolders]):
                                # set status to 1
                                df_ext.loc[df_ext["cmp_id"].astype(int) == int(cmp_id), status_dict[alp_md]] = 1
                            else:
                                # set status to 2
                                df_ext.loc[df_ext["cmp_id"].astype(int) == int(cmp_id), status_dict[alp_md]] = 0

                            #df_ext[status_dict[alp_md]].loc(df_ext["cmp_id"].astype(int) == int(cmp_id)) = 111
            # save df_ext to excel file
            df_ext.to_excel(file, index=False)


    run_alpinac = False
    time_out_alpinac = 30*60 # in s, stops after 15mins / use several hours, if necessary
    run_parallel = True
    if run_alpinac_calc == True:
        if run_parallel == True:
        # Run alpinac for all excel files
            for file in all_excel_files_extended:
                # Read the alpinac input file name
                df_ext = pd.read_excel(Path(file_excel), sheet_name="Sheet1", index_col=0)
                df_alp = df_ext.copy()
                status_list = df_alp["alpinac_status"].to_list()

                cmp_ids_not_yet_processed = df_alp.loc[df_alp["alpinac_status"] == 0, "cmp_id"].to_numpy()
                cmp_ids_extracted_mod = df_alp.loc[(df_alp["extraction_status"] == 1), "cmp_id"].to_numpy()
                cmp_ids = np.intersect1d(cmp_ids_extracted_mod, cmp_ids_not_yet_processed)

                # Initialize alpinac modes
                alpinac_mode_list = ["alpinac_EI_CI", "alpinac_EI_CI_no_N"]

                # dictionaries for alpinac results
                path_for_alpinac_output_dict = {
                    "alpinac_EI_CI": Path(file_excel).parent / "alpinac_results",
                    "alpinac_EI_CI_no_N": Path(file_excel).parent / "alpinac_results_no_N"
                }
                target_elements_dict = {"alpinac_EI_CI": None, "alpinac_EI_CI_no_N": "CHOBrClFIPS"}

                # Create folders for alpinac OUTPUT files if they don't exist
                for path in path_for_alpinac_output_dict:
                    if not path_for_alpinac_output_dict[path].exists():
                        path_for_alpinac_output_dict[path].mkdir()

                # generate alpinac_input_file_mod_dict:
                alpinac_input_file_mod_dict = {}
                for cmp_id in cmp_ids:
                    # read alpinac input file name
                    alp_input_file_relative = df_alp.loc[df_alp["cmp_id"] == cmp_id, "alpinac_input_file_EI_CI"].values[0]
                    # create modified input file
                    alpinac_input_file_mod_dict[cmp_id] = Path(file_excel).parent / "alpinac_input_files" / Path(alp_input_file_relative.replace(".txt", "_mod.txt"))
                # for alpinac_mode in alpinac_mode_list:
                #     # Create a pool of worker processes
                #     with multiprocessing.Pool() as pool:
                #         # Create a shared df
                #         manager = multiprocessing.Manager()
                #         shared_df = manager.list(status_list)  # Create a shared list using Manager

                #         # Define the partial function with shared_df as an argument
                #         partial_process_cmp = partial(process_cmp, shared_df=shared_df)

                #         # Map the partial function to the cmp_ids using multiple processes with a timeout of 20 minutes
                #         results = pool.starmap_async(partial_process_cmp, [(alpinac_input_file_mod_dict[cmp_id], cmp_id, path_for_alpinac_output_dict[alpinac_mode], target_elements_dict[alpinac_mode]) for cmp_id in cmp_ids])
                #         try:
                #             # Wait for the results or timeout after 20 minutes
                #             results.get(timeout=time_out_alpinac)
                #         except multiprocessing.TimeoutError:
                #             # If the timeout is reached, terminate the pool
                #             pool.terminate()
                #         else:
                #             # If the processing completes within the timeout, close and join the pool
                #             pool.close()
                #         pool.join()

                #     # Save the final results
                #     df_alp[status_dict[alpinac_mode]] = shared_df
                #     df_alp.to_excel(file[:-5] + "_extended.xlsx", index=False)

                for alpinac_mode in alpinac_mode_list:
                    #shared_dict = multiprocessing.Manager().dict()  # Create a shared dictionary using Manager

                    # Create a pool of worker processes
                    # Define the partial function with shared_dict as an argument

                    # Use threading to enforce the 20-minute timeout for each process
                    try:
                        result = parallel_alpinac(file_excel)
                        print((f"Process file {file_excel} excecuted."))
                    except:
                        print(f"Process file {file_excel} failed.")
                    #else:
                        #if result is not None:
                        #    shared_dict[cmp_id] = result

                    # Save the final results from the shared_dict to the DataFrame
                    #for cmp_id, result in shared_dict.items():
                    #    idx = np.where(df_alp["cmp_id"].to_numpy() == cmp_id)[0][0]
                    #    df_alp.at[idx, status_dict[alpinac_mode]] = result

                    # df_alp.to_excel(file_excel[:-5] + "_extended.xlsx", index=False)
        else:
            # Run alpinac sequentially for all excel files
            for file in all_excel_files_extended:
                # Read the alpinac input file name
                df_ext = pd.read_excel(Path(file), sheet_name="Sheet1", index_col=0)
                df_alp = df_ext.copy()
                status_list = df_alp["alpinac_status"].to_list()

                cmp_ids_not_yet_processed = df_alp.loc[df_alp["alpinac_status"] == 0, "cmp_id"].to_numpy()
                cmp_ids_extracted_mod = df_alp.loc[(df_alp["extraction_status"] == 1), "cmp_id"].to_numpy()
                cmp_ids = np.intersect1d(cmp_ids_extracted_mod, cmp_ids_not_yet_processed)

                # Initialize alpinac modes
                alpinac_mode_list = ["alpinac_EI_CI", "alpinac_EI_CI_no_N"]

                # dictionaries for alpinac results
                path_for_alpinac_output_dict = {
                    "alpinac_EI_CI": Path(file).parent / "alpinac_results",
                    "alpinac_EI_CI_no_N": Path(file).parent / "alpinac_results_no_N"
                }
                target_elements_dict = {"alpinac_EI_CI": None, "alpinac_EI_CI_no_N": "CHOBrClFIPS"}

                # Create folders for alpinac OUTPUT files if they don't exist
                for path in path_for_alpinac_output_dict:
                    if not path_for_alpinac_output_dict[path].exists():
                        path_for_alpinac_output_dict[path].mkdir()

                # generate alpinac_input_file_mod_dict:
                alpinac_input_file_mod_dict = {}
                for cmp_id in cmp_ids:
                    # read alpinac input file name
                    alp_input_file_relative = df_alp.loc[df_alp["cmp_id"] == cmp_id, "alpinac_input_file_EI_CI"].values[0]
                    # create modified input file
                    alpinac_input_file_mod_dict[cmp_id] = Path(file).parent / "alpinac_input_files" / Path(alp_input_file_relative.replace(".txt", "_mod.txt"))

                # Sequentially process each cmp_id with a time limit of 15 minutes
                for alpinac_mode in alpinac_mode_list:
                    #shared_df = status_list.copy()  # Create a copy of the status_list to simulate shared_df

                    for cmp_id in cmp_ids:
                        # Create a process for each cmp_id
                        process = multiprocessing.Process(target=process_cmp, args=(alpinac_input_file_mod_dict[cmp_id], cmp_id, path_for_alpinac_output_dict[alpinac_mode], target_elements_dict[alpinac_mode]))
                        process.start()
                        
                        # Wait for the process to complete or timeout after 15 minutes
                        process.join(timeout=15 * 60)
                        

                        # If the process is still running after the timeout, terminate it
                        if process.is_alive():
                            process.terminate()
                            process.join()

                    # Save the final results
                    #df_alp[status_dict[alpinac_mode]] = shared_df
                    #df_alp.to_excel(file[:-5] + "_extended.xlsx", index=False)

    update_excel_file_alp_res = True
    if update_excel_file_alp_res == True:
        for file_ext in all_excel_files_extended:
            # check for which cmp_ids the extraction_status is 1
            df_ext = pd.read_excel(Path(file_ext), sheet_name="Sheet1")
            alp_input_filenames = [df_ext.loc[df_ext["cmp_id"] == cmp_id, "alpinac_input_file_EI_CI"].values[0] for cmp_id in df_ext["cmp_id"]]
            # Get input file name
            # Initialize alpinac modes
            alpinac_mode_list = ["alpinac_EI_CI", "alpinac_EI_CI_no_N"]

            # dictionaries for alpinac results
            path_for_alpinac_output_dict = {"alpinac_EI_CI": Path(file_ext).parent / "alpinac_results", "alpinac_EI_CI_no_N": Path(file_ext).parent / "alpinac_results_no_N"}
            target_elements_dict = {"alpinac_EI_CI": None, "alpinac_EI_CI_no_N": "CHOBrClFIPS"}

            for alpinac_mode in alpinac_mode_list:
                # Run alpinac
                path_for_alpinac_output_dir = Path(str(path_for_alpinac_output_dict[alpinac_mode])) #/ (Path(alpinac_input_file_mod.name).stem).replace("input", "output")))
                
                for ionization_mode in ["EI", "CI"]:
                        # get all cmp_ids
                        cmp_ids = df_ext["cmp_id"].to_numpy()
                        # get alpinac results output folders as a dictionary
                        for cmp_id in cmp_ids:
                            # read alpinac input file name
                            alp_input_file_relative = df_ext.loc[df_ext["cmp_id"] == cmp_id, "alpinac_input_file_EI_CI"].values[0]  
                            path_for_alpinac_output_dir_cmp = path_for_alpinac_output_dir / (Path(alp_input_file_relative).stem).replace("input", "output")
                                # get alpinac results as a string
                            flag = False
                            try:
                                string_alpinac_results = get_alpinac_results(path_for_alpinac_output_dir=Path(path_for_alpinac_output_dir_cmp), ionization_mode=ionization_mode)
                                flag = True
                            except:
                                flag = False

                            if flag == False:
                                try:# should have one of the formats, dependent if input file was modified or not
                                    string_alpinac_results = get_alpinac_results(path_for_alpinac_output_dir=Path(str(path_for_alpinac_output_dir_cmp)+"_mod"), ionization_mode=ionization_mode)
                                    flag = True 
                                except:
                                    flag = False
                            
                            if flag == True:
                                print(string_alpinac_results)
                                print("Alpinac ran successfully")            
                                # Set alpinac_status to 1 for cmp_id
                                res_col_name = alpinac_mode.replace("EI_CI", ionization_mode)
                                res_status_col_name = alpinac_mode.replace("_EI_CI", "")  + "_status"
                                df_ext.loc[df_ext["cmp_id"] == cmp_id, res_col_name] = string_alpinac_results
                                df_ext.loc[df_ext["cmp_id"] == cmp_id, res_status_col_name] = 1
                            else:
                                #print("Alpinac failed")
                                # Set alpinac_status to 2 for cmp_id
                                df_ext.loc[df_ext["cmp_id"] == cmp_id, status_dict[alpinac_mode]] = 0
            df_ext.to_excel(file_ext, index=False)

    calculate_match_score_of_alpinac_results = True
    if calculate_match_score_of_alpinac_results == True:
        for file_ext in all_excel_files_extended:
                # check for which cmp_ids the extraction_status is 1
                df_ext = pd.read_excel(Path(file_ext), sheet_name="Sheet1")
                alp_input_filenames = [df_ext.loc[df_ext["cmp_id"] == cmp_id, "alpinac_input_file_EI_CI"].values[0] for cmp_id in df_ext["cmp_id"]]
                # Get input file name
                # Initialize alpinac modes
                alpinac_mode_list = ["alpinac_EI_CI", "alpinac_EI_CI_no_N"]

                # dictionaries for alpinac results
                path_for_alpinac_output_dict = {"alpinac_EI_CI": Path(file_ext).parent / "alpinac_results", "alpinac_EI_CI_no_N": Path(file_ext).parent / "alpinac_results_no_N"}
                target_elements_dict = {"alpinac_EI_CI": None, "alpinac_EI_CI_no_N": "CHOBrClFIPS"}

                for alpinac_mode in alpinac_mode_list:
                    # Run alpinac
                    path_for_alpinac_output_dir = Path(str(path_for_alpinac_output_dict[alpinac_mode])) #/ (Path(alpinac_input_file_mod.name).stem).replace("input", "output")))
                    
                    for ionization_mode in ["EI"]:
                            # get all cmp_ids
                            cmp_ids = df_ext["cmp_id"].to_numpy()
                            # get alpinac results output folders as a dictionary
                            for cmp_id in cmp_ids:
                                # read alpinac input file name
                                alp_input_file_relative = df_ext.loc[df_ext["cmp_id"] == cmp_id, "alpinac_input_file_EI_CI"].values[0]  
                                path_for_alpinac_output_dir_cmp = path_for_alpinac_output_dir / (Path(alp_input_file_relative).stem).replace("input", "output")
                                    # get alpinac results as a string
                                flag = False
                                try:
                                    dict_cs_by_cmp, dict_no_peaks_by_cmp, dict_intensities_summed = get_matchms_results_of_alpinac_fragments(path_for_alpinac_output_dir = Path(path_for_alpinac_output_dir_cmp), database = dbs, target_list = df_peaks, rt_accept_int = (10, 10), ionization_mode = 'EI')
                                    flag = True
                                except:
                                    flag = False

                                if flag == False:
                                    try:# should have one of the formats, dependent if input file was modified or not
                                        dict_cs_by_cmp, dict_no_peaks_by_cmp, dict_intensities_summed = get_matchms_results_of_alpinac_fragments(path_for_alpinac_output_dir=Path(str(path_for_alpinac_output_dir_cmp)+"_mod"), database = dbs, target_list = df_peaks, rt_accept_int = (10, 10), ionization_mode = 'EI')
                                        flag = True 
                                    except:
                                        flag = False

                                #convert dict of dicts to string
                                string_alpinac_cs_results = ""
                                for key in dict_cs_by_cmp:
                                    string_alpinac_cs_results += key + ": " + str(dict_cs_by_cmp[key]) + "\n"        

                                #convert back to dict of dicts -> TODO: in Campaign plotting to extract info! Also import quick_extract functions to generate formula info, e.g. by creatig a new excel
                                dict_cs_by_cmp = {}
                                for line in string_alpinac_cs_results.split("\n"):
                                    if line != "":
                                        key= line.split(": ")[0]
                                        value_str = ":".join(line.split(": ")[1:])
                                        #convert to dict
                                        dict_cs_by_cmp[key] = ast.literal_eval(value_str)                        
                                
                                if flag == True:
                                    print(dict_cs_by_cmp)
                                    print("post alpinac matching ran successfully")            
                                    # Set alpinac_status to 1 for cmp_id
                                    res_col_name = alpinac_mode.replace("EI_CI", "frag_cs")
                                    # declare component with most peaks and their intensities
                                    df_ext.loc[df_ext["cmp_id"] == cmp_id, res_col_name] = string_alpinac_cs_results

                                    dict_main_minor = {0: "main", 1: "minor"}

                                    if len(dict_cs_by_cmp) > 0:
                                        for cmp_ind in range(min(len(dict_cs_by_cmp),2)):
                                            # largest cmp: get first key:
                                            cmp_i = list(dict_no_peaks_by_cmp.keys())[cmp_ind]
                                            res_col_name_largest_comp = alpinac_mode.replace("alpinac_EI_CI", dict_main_minor[cmp_ind] + "_" + str("cs_frag_cmp"))
                                            df_ext.loc[df_ext["cmp_id"] == cmp_id, res_col_name_largest_comp] = str(cmp_i)
                                            res_col_name_largest_comp = alpinac_mode.replace("alpinac_EI_CI", dict_main_minor[cmp_ind] + "_" + str("cs_frag_intens_sum" ))
                                            df_ext.loc[df_ext["cmp_id"] == cmp_id, res_col_name_largest_comp] = str(dict_intensities_summed[cmp_i])
                                            # give largest 0_score and 0_matching_cmp, 1_score and 1_matching_cmp, 2_score and 2_matching_cmp
                                            # 1) get dict for first key
                                            res_score_dict = dict_cs_by_cmp[cmp_i]
                                            for ind, key in enumerate(res_score_dict.keys()):
                                                res_col_name_largest_comp = alpinac_mode.replace("alpinac_EI_CI", dict_main_minor[cmp_ind] + "_" + str("cs_frag_"+str(ind)+"_score"))
                                                df_ext.loc[df_ext["cmp_id"] == cmp_id, res_col_name_largest_comp] = str(res_score_dict[key])
                                                res_col_name_largest_comp = alpinac_mode.replace("alpinac_EI_CI", dict_main_minor[cmp_ind] + "_" + str("cs_frag_"+str(ind)+"_matching_cmp"))
                                                df_ext.loc[df_ext["cmp_id"] == cmp_id, res_col_name_largest_comp] = str(key)
                                                # additional info: formula, iupac_name, cid, inchi, RT_expected
                                                formula_i, iupac_name_i, cid_i, inchi_i, RT_expected_i, formula_dict = formula_by_name(name=clean_empa_name_string(str(key)), formula_dict=formula_dict, formulas_dict_manually_added=formulas_dict_manually_added, target_list=df_peaks_extended)
                                                res_col_name_largest_comp_formula = alpinac_mode.replace("alpinac_EI_CI", dict_main_minor[cmp_ind] + "_" + str("cs_frag_"+str(ind)+"_matching_cmp_formula"))
                                                df_ext.loc[df_ext["cmp_id"] == cmp_id, res_col_name_largest_comp_formula] = str(formula_i)
                                                res_col_name_largest_comp_iupac_name = alpinac_mode.replace("alpinac_EI_CI", dict_main_minor[cmp_ind] + "_" + str("cs_frag_"+str(ind)+"_matching_cmp_iupac_name"))
                                                df_ext.loc[df_ext["cmp_id"] == cmp_id, res_col_name_largest_comp_iupac_name] = str(iupac_name_i)
                                                res_col_name_largest_comp_cid = alpinac_mode.replace("alpinac_EI_CI", dict_main_minor[cmp_ind] + "_" + str("cs_frag_"+str(ind)+"_matching_cmp_cid"))
                                                df_ext.loc[df_ext["cmp_id"] == cmp_id, res_col_name_largest_comp_cid] = str(cid_i)
                                                res_col_name_largest_comp_inchi = alpinac_mode.replace("alpinac_EI_CI", dict_main_minor[cmp_ind] + "_" + str("cs_frag_"+str(ind)+"_matching_cmp_inchi"))
                                                df_ext.loc[df_ext["cmp_id"] == cmp_id, res_col_name_largest_comp_inchi] = str(inchi_i)
                                                res_col_name_largest_comp_RT_expected = alpinac_mode.replace("alpinac_EI_CI", dict_main_minor[cmp_ind] + "_" + str("cs_frag_"+str(ind)+"_matching_cmp_RT_expected"))
                                                df_ext.loc[df_ext["cmp_id"] == cmp_id, res_col_name_largest_comp_RT_expected] = str(RT_expected_i) 

                                    #res_col_name_largest_comp = alpinac_mode.replace("alpinac_EI_CI", "A_frag_cmp_0_score")
                                    #df_ext.loc[df_ext["cmp_id"] == cmp_id, res_col_name_largest_comp] = str(0)


                                    #df_ext.loc[df_ext["cmp_id"] == cmp_id, res_status_col_name] = 1
                                else:
                                    print("matching failed")
                                    # Set alpinac_status to 2 for cmp_id
                                    #df_ext.loc[df_ext["cmp_id"] == cmp_id, status_dict[alpinac_mode]] = 0
                df_ext.to_excel(file_ext, index=False)





        #convert dictionary to columns in dataframe
        
 