## % How to proceed
# 1. do mass calibration for EI and CI using calibration_ei_ci.py
# 2. detect peaks, extract the and compare with Ei data with library entries, write results into excel file
# 3. use interesting compounds, (no carbohydrates, less then 90% confidelity interval) and run alpinac on them
# 4. write results into excel file

#%% Load excel file with library entries
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time


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

logging.basicConfig(level=logging.INFO)
logging.getLogger('alpinac').setLevel(logging.INFO)

def extract_masses_from_mgf(file):
    """Extract masses from mgf file"""
    with open(file, "r") as f:
        for line in f:
            # if start of new spectrum, reset masses
            if line.startswith("BEGIN IONS"):
                masses = []
            # if end of spectrum, return masses
            elif line.startswith("END IONS"):
                return masses
        # otherwise, append mass to masses
            else:
                masses.append(float(line.split()[0].split()[0]))

#%% check if xls file in a subfolder of the current folder
folder_path = r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303"
folder_path = r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303_bad_mass_calib_to_compare"

subfolder_name = "subfolder_name"
extension = "2*results.xlsx"
ionization_for_input = "EI" # uses only EI data, because alpinac is not running well on CI data, TODO change in a later version
run_alpinac = False # will (re)run alpinac, new results will be calculated, even if already available
run_alpinac_EI_CI = False
run_alpinac_no_N = False
use_alpinac = True # alpinac results will be used for results excel file, if available
# CHECK IF DETECT PEAKS WAS ALREADY APPLIED, ELSE (TODO) RUN DETECT PEAKS


# PARAMTERS
cossinsim_limit_for_plotting = 0.9 # assumes that all compounds with a match score larger than this value are identified and plots them as such
rt_range_width = 4 # in seconds, half of the interval around the elution time, where peaks are extracted for alpinac mode extraction
prefiltering = True # prefiltering to run alpinac? True or False
min_peak_no_to_run_alpinac = 3 # minimum number of peaks extracted to run alpinac
match_score_to_run_alpinac = [0,1] # run alpinac on all compounds with a match score between a and b (all) b=1 : 100%, b=0 0%
plot_summary_pie = True # plot summary pie chart? True or False

# CHECK IF FOLDER BY DETECT PEAKS WAS CREATED, AND LIST ALL EXCEL FILES

# get a list of all files and folders in the specified directory
dir_contents = os.listdir(folder_path)

# iterate over each item in the directory
for item in dir_contents:
    # construct the full path to the item
    item_path = os.path.join(folder_path, item)
    
    # check if the item is a folder
    if os.path.isdir(item_path):
        print("processing folder: " + item)
        print(item_path)
        # construct the full search path using the folder path and extension
        search_path = os.path.join(item_path, extension)

        # use the glob module to find all files that match the search path
        excel_files = glob.glob(search_path)

        # print the list of excel files
        print(excel_files)

#%%
print(excel_files)

#%% OPEN EXCEL FILE AND READ OUT POTENTIALLY INTERESTING COMPOUNDS
from pathlib import Path

for file in excel_files:
    print("processing file: " + file)
    # try open excel file
    try:
        df = pd.read_excel(Path(file), sheet_name="Sheet1")
        print("excel file opened")
        # add info of results of extraction with mygu script
        # check if copy of excel file exists
        if os.path.isfile(file[:-5] + "_mod_kaho.xlsx"):
            print("modifed excel file already exists")
            df_copy = pd.read_excel(Path(file[:-5] + "_mod_kaho.xlsx"), sheet_name="Sheet1")
        else:
            print("modifed excel file does not exist, create new one")
            # make copy of excel file
            df_copy = df.copy()
            # add empty column
            df_copy["Extraction"] = ""
            # read score values of column "0_score" into a numpy array
            df_copy["alpinac_results"] = ""
            # read score values of column "0_score" into a numpy array
            df_copy["alpinac_results_no_N"] = ""
            # read score values of column "0_score" into a numpy array
            df_copy["alpinac_results_EI_CI"] = ""
            # save as new excel file
            df_copy.to_excel(file[:-5] + "_mod_kaho.xlsx", index=False)


        scores = df["0_score"].to_numpy()
        # read prominences of column "total_prominences" into a numpy array
        prominences = df["total_prominences"].to_numpy()
        # number of entries smaller than cossinsim_limit_for_running_alpinac (unidentified)
        print("number of compounds smaller than cossinsim_limit_for_running_alpinac: " + str(np.sum(scores > cossinsim_limit_for_plotting)))
        ind_smaller_09 = np.where(scores < cossinsim_limit_for_plotting)
        # number of entries larger than cossinsim_limit_for_running_alpinac 
        print("number of compounds larger than cossinsim_limit_for_running_alpinac: " + str(np.sum(scores > cossinsim_limit_for_plotting)))
        # index of entries larger than cossinsim_limit_for_running_alpinac 
        ind_larger_09 = np.where(scores >= cossinsim_limit_for_plotting)
        # print compounds in column df["0_matching_cmp"] with value (cosine similarity) larger than cossinsim_limit_for_running_alpinac
        print("compounds with cosine similarity larger than cossinsim_limit_for_running_alpinac, ordered by prominence: " + str(np.array(df["0_matching_cmp"])[ind_larger_09]))
        # index of entries larger than cossinsim_limit_for_running_alpinac and prominence larger than 0.1 of sum of all prominences
        ind_larger_09_prom = np.where((scores >= cossinsim_limit_for_plotting) & (prominences > 0.02*np.sum(prominences)))
        import matplotlib.pyplot as plt
        # generate some sample data
        prominences = np.array(df["total_prominences"])[ind_larger_09_prom]
        # sum of < 2% identified compounds
        prominences = np.append(prominences, np.sum(np.array(df["total_prominences"])[ind_larger_09]) -  np.sum(np.array(df["total_prominences"])[ind_larger_09_prom]))
        prominences = np.append(prominences, np.sum(np.array(df["total_prominences"])[ind_smaller_09]))

        labels_tot = np.array(df["0_matching_cmp"])[ind_larger_09_prom]
        # use the first part before "_" of labels for the pie chart
        labels = [label.split("_")[0] for label in labels_tot]
        labels.append("identified < 2%")
        labels.append("unidentified")
        # the second part of label tot use for the origin (NIST, mygu)
        origin = [label.split("_")[1] for label in labels_tot]
        origin.append("identified")
        origin.append("unidentified")

        if plot_summary_pie == True:
            #plot two pie charts next to each other
            fig, (ax1, ax2) = plt.subplots(2, 1)
            # create a figure and axis for the plot
            #fig, ax = plt.subplots()
            # create the pie chart
            wedges, labels, _ = ax1.pie(prominences, labels=labels, autopct='%1.1f%%', labeldistance=1.1, pctdistance=0.85, radius = 1,  textprops={'fontsize': 8})
            # textcolor is dark blue if origin is Myriam, else black
            for label_ind in range(len(labels)):
                if origin[label_ind] == "myriam":
                    labels[label_ind].set_color('blue')
            # ast, unknown wedge is white
            wedges[-1].set_color('lightgrey')
            # add name of the file as title
            title_pie = str(file.split("\\")[-1].split("_")[0])
            ax1.set_title(title_pie)
            #ax.set_title("Prominences by Compound")
            plt.tight_layout()
            # show the plot
            #plt.show()

            #pie plot of < 2% identified compounds
            ind_larger_09_small_prom = np.where(np.array(df["0_score"].to_numpy() >= cossinsim_limit_for_plotting) & (np.array(df["total_prominences"]) <= 0.02*np.sum(np.array(df["total_prominences"]))))
            prominences_small_prom = np.array(df["total_prominences"])[ind_larger_09_small_prom]
            labels_tot = np.array(df["0_matching_cmp"])[ind_larger_09_small_prom]
            # use the first part before "_" of labels for the pie chart
            labels_small_prom = [label.split("_")[0] for label in labels_tot]
            # the second part of label tot use for the origin (NIST, mygu)
            origin_small_prom = [label.split("_")[1] for label in labels_tot]
            # generate greyshade color array for the pie chart
            darkgrey_colors = [plt.cm.Greys(ind) for ind in np.linspace(0.5, 0.75, len(labels_small_prom))]

            # create the pie chart
            wedges, labels_small_prom, _ = ax2.pie(prominences_small_prom[0:30], labels=labels_small_prom[0:30], autopct='%1.1f%%', colors = darkgrey_colors, labeldistance=1.1, pctdistance=0.85, radius = 1,  textprops={'fontsize': 8})
            # textcolor is dark blue if origin is Myriam, else black
            for label_ind in range(len(labels_small_prom)):
                if origin_small_prom[label_ind] == "myriam":
                    labels_small_prom[label_ind].set_color('blue')
            # add name of the file as title
            ax2.set_title("Rel. prominences of identified <2%")
            plt.tight_layout()
            # show the plot
            # get the text objects from the pie chart
            #texts = plt.gca().texts

            # set the fontsize of the text objects
            #for text in texts:
            #    text.set_fontsize(6)
            plt.show()

            # save pie plot
            fig.savefig(file.split(".xlsx")[0] + "_pie_overview.png", dpi=150)


            # Show contribution of largest unidentified compounds
            fig, (ax1) = plt.subplots(1)
            #pie plot of < 2% identified compounds
            prominences_unidentified = np.array(df["total_prominences"])[ind_smaller_09]
            labels_unidentified = [str(str(ind)+":ID "+str(df["cmp_id"][ind])+", RT="+str(round(df["RT"][ind]))+", "+str(df["number_of_peaks"][ind])+" peaks") for ind in ind_smaller_09[0]]
            # alert for too little peaks
            origin_unidentified = [df["number_of_peaks"][ind] for ind in ind_smaller_09[0]]
            # generate greyshade color array for the pie chart
            darkgrey_colors = [plt.cm.Greys(ind) for ind in np.linspace(0.1, 0.5, len(labels_unidentified))]
            # create the pie chart
            wedges, labels_unidentified, _ = ax1.pie(prominences_unidentified[0:30], labels=labels_unidentified[0:30], autopct='%1.1f%%', colors = darkgrey_colors, labeldistance=1.1, pctdistance=0.85, radius = 1,  textprops={'fontsize': 8})
            # textcolor is dark blue if origin is Myriam, else black
            for label_ind in range(len(labels_unidentified)):
                if origin_unidentified[label_ind] <=3:
                    labels_unidentified[label_ind].set_color('red')
            # add name of the file as title
            ax1.set_title("Rel. prominences of unidentified")
            plt.tight_layout()
            # show the plot
            plt.show()

            fig.savefig(file.split(".xlsx")[0] + "_pie_unidentified.png", dpi=150)

        # PREFILTERING

        if prefiltering == True:
            min_peak_no_to_run_alpinac = 3 # lim = 0
            match_score_to_run_alpinac = [0,1] # run alpinac on all compounds with a match score between 0 and 1 (all)
            ind_sel_by_score_lim = np.intersect1d(np.where(match_score_to_run_alpinac[0] <= scores), np.where(scores <= match_score_to_run_alpinac[1]))
            ind_sel_by_min_no_pks = np.where(np.array(df["number_of_peaks"]) >= min_peak_no_to_run_alpinac)
            # indeci where both conditions are fullfilled
            ind_sel = np.intersect1d(ind_sel_by_score_lim, ind_sel_by_min_no_pks)
        else:
            # all indici
            ind_sel = np.arange(len(df["number_of_peaks"]))
        # Print the compounds that are going to be processed
        print("compounds that are going to be processed: " + str(np.array(df["cmp_id"])[ind_sel]))


        # Prefilter candidates which are not identified, use only the ind_smaller_09 elements for which the number of peaks is larger than 3
        #ind_smaller_09_sel = ind_smaller_09[0][np.array(df["number_of_peaks"])[ind_smaller_09] > 3]
        # get the mgf files from the folder matchms
        matchms_dir = Path(file).parent /"matchms_spectra"
        # create the corresponding list of matchms files
        matchms_files = [matchms_dir/str(str(df["cmp_id"][ind])+'.mgf') for ind in ind_sel]



        # try using alpinac mode extraction to extract the peaks at given retention times for the compounds of interest
   

        # select indeci
        indici_process = ind_sel#[0:2]
        for ind_proc in indici_process:
            print("processing compound " + str(df["cmp_id"][ind_proc]) + " at elution time " + str(df["RT"][ind_proc]))

            # read excelfile into df_copy, which will be modified by the script
            df_copy = pd.read_excel(Path(file[:-5] + "_mod_kaho.xlsx"), sheet_name="Sheet1")

            # get the retention times of the compounds of interest
            elution_time = np.array(df["RT"])[ind_proc]
            # limit to masses found with python peak finder (better performance)
            # going to corresponding matchms file:
            # try to get the matchms files:
            try:
                # get the mgf files from the folder matchms
                matchms_dir = Path(file).parent /"matchms_spectra"
            except:
                print("Error in matchms file path, check if matchms folder was copied to measurement folder")
                #continue
            # create the corresponding list of matchms files
            matchms_file = matchms_dir/str(str(df["cmp_id"][ind_proc])+'.mgf')

            # extract the masses from the mgf files
            masses = extract_masses_from_mgf(matchms_file)
            # convert masses entries to integer, because Myriams script uses integers
            masses = [int(mass) for mass in masses]

            raw_data_file = Path(file).parent.parent /str(Path(file).parent.name + ".h5")

            #

            #dict_spikes = {}

            #file = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230401.1747.air.3.h5")

            #segments=(6,10),
            #segments=(1,2),
            ionization_segments={
                IonizationMethod.EI: (3,4),
                IonizationMethod.CI: (1,2),
            }


            #i=2
            #cmd_options = ["--rt-start", "--rt-stop", "--mass-start", "--mass-stop", "--mass-list"]
            #elution_time = 1654
            

            alpinac_input_file = Path(file).parent /str(Path(file).parent.name + '.frag.'+str(int(elution_time -rt_range_width))+'s'+ str(int(elution_time +rt_range_width))+'s'+'.txt')
            #alpinac_input_file = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230401.1747.air.3/230401.1747.air.3.frag.1647s1661s_EI_only_min_4_peaks_per_compound_mod.txt")
            # GENERATE ALPINAC INPUT FILE
            try:
                if alpinac_input_file.exists():
                    print("input file already exists, skipping")
                    #continue
                else:
                    print("alpinac input file does not exist, creating")
                    nontarget_peak_extraction(
                        path_file = raw_data_file,                              
                        rt_start_extract = elution_time-rt_range_width,                              
                        rt_stop_extract = elution_time+rt_range_width, 
                        adduct= "H", 
                        #plot_pseudovoigt = True,
                        #coelution_plot_mode = True,
                        #mass_bin_plot_mode = True,
                        #save_plots=[Plots.CHROMATOGRAPHY_FITTING, Plots.COELUTION_PLOT, Plots.RT_PEAK_DETECTION],
                        #mass_start = mass_start,                              
                        #mass_stop = mass_stop,                              
                        #mass_list = masses,
                        ionization_segments=ionization_segments,
                        rt_sigma_guess_scale=1.0, # might give wrong result for large RT 
                        mass_u_ppm_min = 10.0,
                    )
                    print("input file generated")
            except:
                print("Error in alpinac mode extraction")
                #continue
            
            #analyze the result file
            # get the name of the input file


            # RUN ALPINAC
            try:
                # read input file
                df_input = pd.read_csv(alpinac_input_file, delim_whitespace=True)
                #get ionization mode
                ionization_mode = df_input["Ionisation"]
                # if there are no EI peaks, skip the file
                if len(df_input[df_input["Ionisation"] == "EI"]) == 0:
                    print("no EI peaks found, skipping")
                    counts_per_compound_text = "no EI peaks extr."
                else:
                    print("modifing alpinac input file")

                    # modify the input file
                    # delete CI data, TODO change in a later version, so that CI data is also used
                    if ionization_for_input == "EI":
                        df_input = df_input[df_input["Ionisation"] == "EI"]    
                    # if to litte EI peaks for one compound, delete the compound
                    # loop over all non-negative compound bins and keep only the ones with more than 3 peaks
                    for i in np.unique(df_input["compound_bin"]):
                        # if there are less than 3 peaks, delete the compound
                        if len(df_input[df_input["compound_bin"] == i]) < 3:
                            df_input = df_input[df_input["compound_bin"] != i]
                    # convert all NaN values to "None"
                    df_input = df_input.fillna("None")
                    # make sure that Adduct column is a string
                    #df_input["Adduct"] = df_input["Adduct"].astype(str)

                    # save the modified file
                    modified_alpinac_input_file_path = alpinac_input_file.parent /(alpinac_input_file.stem +"_EI_only"+alpinac_input_file.suffix)
                    #save as text file

                    df_input.to_csv(modified_alpinac_input_file_path, sep="\t", index=True, header=True)
                    # print how many occurences of each compound bin were found
                    counts_per_compound = df_input["compound_bin"].value_counts()
                    # remove counts belongign to negative compound bins
                    counts_per_compound = counts_per_compound[counts_per_compound.index >= 0]
                    if len(counts_per_compound) == 0:
                        counts_per_compound_text = "no good EI peaks extr."
                    else:
                        counts_per_compound_text =  [str(item) +"("+str(counts_per_compound[item])+")" for item in counts_per_compound.index.tolist()]
                        # flatten list of strings
                        counts_per_compound_text = "#compound (#peaks): " + str(counts_per_compound_text).replace("[", "").replace("]", "").replace("'", "")
                        # save a file suitable for alpinac onlyconsidering minimum compund with more than 3 peaks
                        df_input_suitable = df_input[df_input["compound_bin"].isin(counts_per_compound[counts_per_compound > 3].index.tolist())]
                        modified_alpinac_input_file_path = alpinac_input_file.parent /(alpinac_input_file.stem +"_EI_only_min_4_peaks_per_compound"+alpinac_input_file.suffix)
                        df_input_suitable.to_csv(modified_alpinac_input_file_path, sep="\t", index=False)
                    print(counts_per_compound_text )
                    print("alpinac input file modified")




                df_copy["Extraction"][ind_proc] = counts_per_compound_text
                # save the copy
                # open excel file and add one value and close it
    
                

                if use_alpinac:
                    #test if inputfile has data:
                    # read input file
                    #df_input = pd.read_csv(modified_alpinac_input_file_path, delim_whitespace=True)
                    #df_input = pd.read_csv(alpinac_input_file, delim_whitespace=True); modified_alpinac_input_file_path = alpinac_input_file
                    if len(df_input) == 0:
                        print("Not enough EI peaks found, skipping alpinac")
                        #continue
                    else:
                        # py -m alpinac.mode_identification C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230401.1747.air.3\230401.1747.air.3.frag.1632s1652s_EI_only.txt --fitting-method discrete -d -r
                        # run alpinac
                        try:
                            # open excel file
                            for alp_i in range(3):
                                
                                if alp_i == 0:
                                    #path_for_alpinac_output = Path(modified_alpinac_input_file_path.parent) / (modified_alpinac_input_file_path.stem)
                                    path_for_alpinac_output = alpinac_input_file.parent /(alpinac_input_file.stem +"_EI_only_min_4_peaks_per_compound")
                                    if run_alpinac == True:
                                        print("try running alpinac with all atoms")
                                        # run alpinac with all atoms
                                        make_identification_from_file(path_file=modified_alpinac_input_file_path, output_dir = path_for_alpinac_output, show_plt_figures= False) 
                                        print("alpinac run wit all possible atoms")
                                    else:
                                        print("precalulated alpinac results are used, if available")
                                    col_name_results = "alpinac_results"
                                    suffix = ""
                                elif alp_i == 1:
                                    path_for_alpinac_output = Path(modified_alpinac_input_file_path.parent) / (modified_alpinac_input_file_path.stem + "_no_N")
                                    if run_alpinac_no_N == True:
                                        print("try running alpinac with all atoms")
                                        # run alpinac with all atoms but N
                                        # run alpinac with only C,H,N,O,S
                                        make_identification_from_file(path_file=modified_alpinac_input_file_path, output_dir = path_for_alpinac_output, show_plt_figures= False, target_elements="CHOBrClFIPS") 
                                        print("alpinac run without N")
                                    else:
                                        print("precalulated alpinac results are used, if available")
                                    suffix = "_no_N"
                                    col_name_results = "alpinac_results_no_N"
                                else: 
                                    path_for_alpinac_output= Path(alpinac_input_file.parent) / (alpinac_input_file.stem + "_EI_CI")    
                                    # run alpinac including CI data
                                    if run_alpinac_EI_CI == True:
                                        print("try running alpinac with all atoms")
                                        # run alpinac with all atoms but N
                                        # run alpinac with only C,H,N,O,S
                                        make_identification_from_file(path_file=alpinac_input_file, output_dir = path_for_alpinac_output, show_plt_figures= False, target_elements="CHOBrClFIPS") 
                                        print("alpinac run for EI, CI")
                                    else:
                                        print("precalulated alpinac results are used, if available")
                                    suffix = "_no_N"
                                    col_name_results = "alpinac_results_EI_CI"

                                        
                                #continue
                                    
                                # add alpinac results to the excel file
                                # read alpinac results
                                # list all directories in Path
                                #path_for_alpinac_output = Path(modified_alpinac_input_file_path.parent) / Path(modified_alpinac_input_file_path.stem + suffix)
                                dir_contents = os.listdir(path_for_alpinac_output)
                                # dir_contents = os.listdir(path_for_alpinac_output = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230401.1747.air.3\230401.1747.air.3.frag.1828s1848s_EI_only"))
                                # get the directories in path starting with "Compound_"
                                compound_dirs = list(filter(lambda x: x.startswith("Compound_"), dir_contents))
                                string_alpinac_results = ""
                                for compound_dir in compound_dirs:
                                    # open most_likely_molecular_ion.txt in the compound dir as pandas dataframe
                                    df_most_likely_molecular_ion = pd.read_csv(path_for_alpinac_output / compound_dir / "most_likely_mol_ions.txt", delim_whitespace=True)
                                    #add column names
                                    df_most_likely_molecular_ion.columns =  ["spectrum_id", "max_adduct", "ranking", "formula", "DBE", "likelihood", "", ""]
                                    # delete rows if  df_most_likely_molecular_ion["spectrum_id"] != max(df_most_likely_molecular_ion["spectrum_id"])   # The largest one (1) should be EICI the first EI only
                                    #ind_val = np.where(df_most_likely_molecular_ion["spectrum_id"] == np.max(df_most_likely_molecular_ion["spectrum_id"]))
                                    #df_most_likely_molecular_ion = pd.datanp.array(df_most_likely_molecular_ion)[ind_val]
                                    #df_most_likely_molecular_ion = df_most_likely_molecular_ion[np.where(df_most_likely_molecular_ion["spectrum_id"] == max(df_most_likely_molecular_ion["spectrum_id"]))]                                    
                                    # get the most likely molecular ion
                                    #get max compound

                                    string_alpinac_results += compound_dir + ": "
                                    for i in range(np.min([3, len(df_most_likely_molecular_ion)])):
                                        print(i)
                                        if df_most_likely_molecular_ion["spectrum_id"][i] == np.max(df_most_likely_molecular_ion["spectrum_id"]):
                                            most_likely_molecular_ion = df_most_likely_molecular_ion["formula"][i]
                                            most_likely_molecular_ion_score = df_most_likely_molecular_ion["likelihood"][i]
                                            rank = df_most_likely_molecular_ion["ranking"][i]
                                            string_alpinac_results += "rank "+ str(rank) + ": " + str(most_likely_molecular_ion) + " (" + str(most_likely_molecular_ion_score) + "%) "
                                            string_alpinac_results += "; "
                                # add alpinac results to the excel file
                                df_copy[col_name_results][ind_proc] = string_alpinac_results
                                # save the copy
                                # df_copy.to_excel(Path(file).parent /str(Path(file).parent.name + "_mod_kaho.xlsx"), index=False)
                        except:
                            print("Error in alpinac mode extraction")
                            continue


                #save results in excel file

                df_copy.to_excel(Path(file[:-5] + "_mod_kaho.xlsx"), index=False)

            except:
                print("Error in alpinac mode extraction")
                continue

    
    except:
        print("excel file could not be opened, check rights")
        continue



#df = pd.read_excel(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230401.1747.air.3\230401.1747.air.3_peak_identification_results.xlsx", sheet_name="Sheet1")

