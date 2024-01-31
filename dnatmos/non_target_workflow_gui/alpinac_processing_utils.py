from pathlib import Path
import pandas as pd
import numpy as np
from db_reader import AlpinacFragmentsReader, JdxReader, AbstractReader
# for post-matching:
from matchms.Spectrum import Spectrum
from matchms.similarity import CosineGreedy
from matchms.filtering import normalize_intensities
from alpinac_gui.matchms_funcs import formula_from_namestring, formula_to_latex, halocarbname2formula, plot_mass_spectra, load_from_mgf, metadata_processing, peak_processing, AlpinacData, halocarbname2formula




def create_modified_input_file(alp_input_filename: Path, extended_quick_res_path: Path, df_ext: pd.DataFrame, cmp_id: int):
    ind_cmp = df_ext.index[df_ext["cmp_id"] == cmp_id].to_numpy()[0]    
                       # try to open input file
    try:
        input_df = pd.read_csv(alp_input_filename, delim_whitespace=True)
    except:
        print("Error in reading alpinac input file {}".format(alp_input_filename))
        # check if file exists
        if not alp_input_filename.exists():
            print("alpinac input file {} does not exist".format(alp_input_filename))
            return df_ext
    if sum(input_df["compound_bin"]>=0) < 1:
        df_ext.loc[ind_cmp, "extraction_status"] = 2 # set to failed
        print("only noise extracted")
    else:
        cmp_bins_unique = input_df["compound_bin"].unique()
        # remove -1 from unique bins if present
        if -1 in cmp_bins_unique:
            cmp_bins_unique = cmp_bins_unique[cmp_bins_unique != -1]
        if len(cmp_bins_unique) < 1:
            df_ext.loc[ind_cmp, "extraction_status"] = 2 # set to failed
        else:
        # go through all positive ints other in compound_bin and count the occurences
        # if there are more than 3 in at least one, set extraction_status to 1
        # if there are less than 3 in all, set extraction_status to 2
            EI_occurences = []
            CI_occurences = []
            for cmp_bin in cmp_bins_unique:
                #sum occurences for ionization = "EI" and "CI"
                indici_cmp_bin_occurrences_EI = sum((input_df["compound_bin"] == cmp_bin) & (input_df["Ionisation"] == "EI"))
                EI_occurences.append(indici_cmp_bin_occurrences_EI)
                indici_cmp_bin_occurrences_CI = sum((input_df["compound_bin"] == cmp_bin) & (input_df["Ionisation"] == "CI"))
                CI_occurences.append(indici_cmp_bin_occurrences_CI)
            # if there are more than 3 in at least one, set extraction_status to 1
            if (max(EI_occurences) > 3) or (max(CI_occurences) > 3):
                df_ext.loc[ind_cmp, "extraction_status"] = 1 # set to extracted
            else:
                df_ext.loc[ind_cmp, "extraction_status"] = 3
            # set extraction_number_compounds
            df_ext.loc[ind_cmp, "extraction_number_cmps"] = str(cmp_bins_unique)
            df_ext.loc[ind_cmp, "extraction_peaks_per_cmp_EI"] = str(EI_occurences)
            df_ext.loc[ind_cmp, "extraction_peaks_per_cmp_CI"] = str(CI_occurences)
            input_df_mod = input_df.copy()
            # set all compound_bins to -1 if there are less than 3 in EI_occurences
            for cmp_bin, occurence in zip(cmp_bins_unique, EI_occurences):
                if occurence < 3:
                    input_df_mod.loc[input_df_mod["compound_bin"] == cmp_bin, "compound_bin"] = -1
            # save modified input file
            alp_input_filename_mod = Path(str(alp_input_filename)[:-4] + "_mod.txt")
            input_df_mod.where(pd.notnull(input_df_mod), 'None').to_csv(alp_input_filename_mod,  index=False, sep="\t")
        # save df_ext to excel file
        df_ext.to_excel(extended_quick_res_path, index=False)
    return df_ext


def get_alpinac_results(path_for_alpinac_output_dir, ionization_mode = "EI"):                       
    # add alpinac results to the excel file
    dir_contents = os.listdir(path_for_alpinac_output_dir)
    # dir_contents = os.listdir(path_for_alpinac_output = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230401.1747.air.3\230401.1747.air.3.frag.1828s1848s_EI_only"))
    # get the directories in path starting with "Compound_"
    compound_dirs = list(filter(lambda x: x.startswith("Compound_"), dir_contents))
    string_alpinac_results = ""
    for compound_dir in compound_dirs:
        result_most_likely_ions_file = Path(path_for_alpinac_output_dir) / Path(compound_dir) / "most_likely_mol_ions.txt"
        # chcek if file exists
        if result_most_likely_ions_file.exists():
            # open most_likely_molecular_ion.txt in the compound dir as pandas dataframe
            df_most_likely_molecular_ion = pd.read_csv(result_most_likely_ions_file, delim_whitespace=True)
            #add column names
            df_most_likely_molecular_ion.columns =  ["spectrum_id", "max_adduct", "ranking", "formula", "DBE", "likelihood", "", ""]
            # filter by ionization mode
            if ionization_mode == "EI":
                df_most_likely_molecular_ion = df_most_likely_molecular_ion.loc[df_most_likely_molecular_ion["spectrum_id"] == 0]
            elif ionization_mode == "CI":
                df_most_likely_molecular_ion = df_most_likely_molecular_ion.loc[df_most_likely_molecular_ion["spectrum_id"] == 1]
            else:
                raise Exception("ionization_mode must be EI or CI")
            #get max compound
            string_alpinac_results += compound_dir + ": "
            for i in range(np.min([3, len(df_most_likely_molecular_ion)])):
                #print(i)
                #if df_most_likely_molecular_ion["spectrum_id"][i] == np.max(df_most_likely_molecular_ion["spectrum_id"]):
                most_likely_molecular_ion = df_most_likely_molecular_ion["formula"].values[i]
                most_likely_molecular_ion_score = df_most_likely_molecular_ion["likelihood"].values[i]
                rank = df_most_likely_molecular_ion["ranking"].values[i]
                string_alpinac_results += "rank "+ str(rank) + ": " + str(most_likely_molecular_ion) + " (" + str(most_likely_molecular_ion_score) + "%) "
                string_alpinac_results += "; "
    # add alpinac results to the excel file
    return string_alpinac_results


def get_matchms_results_of_alpinac_fragments(path_for_alpinac_output_dir, database, target_list, rt_accept_int=(10, 10), ionization_mode='EI'):
    """
    This function compares the alpinac fragments with the database and returns the results in a dictionary.
    
    Parameters
    ----------
        path_for_alpinac_output_dir : str
        Path to the alpinac output directory.
        database : dict
        Dictionary with the database entries.
        target_list : pandas.DataFrame
        Dataframe with the target list.
        rt_accept_int : tuple
        Tuple with the RT acceptance interval (left, right).
        ionization_mode : str
    """
    cosine_greedy = CosineGreedy(tolerance=0.2)

    # add alpinac results to the excel file
    # path_for_alpinac_output_dir = str(Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303_interesting_file\230311.0125.tank.1\alpinac_results\230311.0125.tank.1_EI_CI_output_198_mod"))
    dir_contents = os.listdir(path_for_alpinac_output_dir)
    # dir_contents = os.listdir(path_for_alpinac_output_dir = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230401.1747.air.3\230401.1747.air.3.frag.1828s1848s_EI_only\230311.0125.tank.1_EI_CI_output_11_mod"))
    # get the directories in path starting with "Compound_"
    compound_dirs = list(filter(lambda x: x.startswith("Compound_"), dir_contents))
    string_alpinac_results = ""
    df_peaks = target_list.copy()
    time_factor_neg = rt_accept_int[0]
    time_factor_pos = rt_accept_int[1]
    dict_res_cs = {}
    dict_no_peaks_cmp = {}
    dict_intensities_summed = {}


    for compound_dir in compound_dirs:
        res_file = Path(path_for_alpinac_output_dir) / Path(compound_dir) / "results_file_mygu.txt"
        # chcek if file exists
        if res_file.exists():
            # get alpinac spectrum from results_file_mygu.txt
            alp_spec = AlpinacData.from_alpinac(res_file, ionization_mode = "EI")
            # zero value will be ignored when matching (in range only)
            # mz
            mz = alp_spec.peaks.mz
            i = alp_spec.peaks.intensities
            # get sum of intensities
            intensities_summed = np.sum(i)
            key = compound_dir.replace("Compound_", "c")
            dict_intensities_summed[key] = intensities_summed

            no_peaks = len(mz)
            key = compound_dir.replace("Compound_", "c")
            dict_no_peaks_cmp[key] = no_peaks


            # get average RT of alpinac spectrum
            rt = np.mean(alp_spec.rt_in_s)
            # get RT window of alpinac spectrum
            cmp_start_rt = np.min(alp_spec.rt_in_s)
            cmp_end_rt = np.max(alp_spec.rt_in_s)


            
            # calculate matching results (see detect_peaks)
            # TODO: put code below in a function
            # Get compounds matching with the rt
            mask_matching_time = (
                cmp_start_rt > df_peaks["RT"] - time_factor_neg * df_peaks["RT-Window"] / 2
            ) & (cmp_end_rt < df_peaks["RT"] + time_factor_pos * df_peaks["RT-Window"] / 2)
            compounds = list(df_peaks.loc[mask_matching_time, "Substance"])

            #matchms_spectra = Spectrum(mz=np.array(mz, dtype=float), intensities=i)


            #with open(spectra_dir / f"{cmp_id}.mgf", 'w') as f:
            #    save_as_mgf([matchms_spectra], f)
            #ms_spectras.append(matchms_spectra)
            # Compare the compound to known compounds

            scores = {}
            matches = {}

            for compound in compounds:

                for db, db_compounds in database.items():

                    if compound not in db_compounds:
                        #print("Unknown RT of", compound, "in db", db)
                        continue
                    # Compare the two spectra
                    reference = db_compounds[compound]
                    # use only masses in the range of the experimental device
                    mask_reference_in_range = (reference.peaks.mz > 27) & (reference.peaks.mz < 500)
                    reference_in_range  = Spectrum(mz=reference.peaks.mz[mask_reference_in_range], intensities=reference.peaks.intensities[mask_reference_in_range], metadata={"precursor_mz": -1})
                    score_res = cosine_greedy.pair(
                        normalize_intensities(reference_in_range), normalize_intensities(alp_spec)
                    )
                    w_cossim = 0.5
                    w_no_peaks = 0.5
                    
                    if compound not in matches:
                        matches[compound] = {}
                    score_cossim = float(score_res["score"])
                    matched_peaks = float(score_res["matches"])
                    # total number of peaks in the reference
                    total_no_peaks = len(reference_in_range.peaks)
                    # get all reference peaks indici with mz = 30 and intensity > 0.05
                    indici = np.argwhere(reference_in_range.peaks.intensities > 0.05)
                    total_no_peaks_in_lim = len(indici)

                    # score is product of cosine similarity (intensity match) and relative number of matched peaks (reliabiltiy of match)
                    if score_cossim == 0:
                        score = 0 # avoid zero division error
                    else:
                        # weighted ratio (by intensity) of matched peaks to total peaks
                        weighted_ratio = matched_peaks/total_no_peaks
                        if total_no_peaks_in_lim < matched_peaks:
                            score = w_cossim*score_cossim 
                            #logger.warning(f"More peaks matched than expected in measured spectra (mz>27, peak hight 5% of max), adapt measured range, in the reference for {compound} in {db}")
                        else:
                            score = w_cossim*score_cossim + w_no_peaks * matched_peaks/total_no_peaks_in_lim
                    # Add score to the results dictionaries
                    scores[f"{compound}_{db}"] = score
                    matches[compound] = score

                    # get largest 3 scores:
                    # sort scores by value
                    sorted_scores = {k: v for k, v in sorted(scores.items(), key=lambda item: item[1], reverse=True)}
                    # get first 3 keys
                    sorted_scores_keys = list(sorted_scores.keys())
                    sorted_scores_keys = sorted_scores_keys[:3]
                    # get first 3 values
                    sorted_scores_values = list(sorted_scores.values())
                    sorted_scores_values = sorted_scores_values[:3]
                    # get first 3 scores
                    sorted_scores = dict(zip(sorted_scores_keys, sorted_scores_values))
                    print(sorted_scores)
                    key = compound_dir.replace("Compound_", "c")
                    dict_res_cs[key] = sorted_scores

                    # order dicts by value of dict_no_peaks
                    dict_no_peaks_cmp = {k: v for k, v in sorted(dict_no_peaks_cmp.items(), key=lambda item: item[1], reverse=True)}
                    # get sorted keys
                    sorted_keys = list(dict_no_peaks_cmp.keys())
                    # order dict_res_cs by sorted keys
                    dict_res_cs = {k: dict_res_cs[k] for k in sorted_keys}
                

    return dict_res_cs, dict_no_peaks_cmp, dict_intensities_summed


if __name__ == "__main__":
    #path_for_alpinac_output_dir = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303_interesting_file\230311.0125.tank.1\alpinac_results\230311.0125.tank.1_EI_CI_output_198_mod")
    #path_for_alpinac_output_dir = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303_interesting_file\230311.0125.tank.1\alpinac_results\230311.0125.tank.1_EI_CI_output_198")
    path_for_alpinac_output_dir = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303_interesting_file\230311.0125.tank.1\alpinac_results\230311.0125.tank.1_EI_CI_output_198_mod")
    #path_for_alpinac_output_dir = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303_interesting_file\230311.0125.tank.1\alpinac_results\230311.0125.tank.1_EI_CI_output_198_mod")
    #path_for_alpinac_output_dir = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303_interesting_file\230311.0125.tank.1\alpinac_results\230311.0125.tank.1_EI_CI_output_198_mod")
    #path_for_alpinac_output_dir = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303_interesting_file\230311.0125.tank.1\alpinac_results\230311.0125.tank.1_EI_CI_output_198_mod")
    #path_for_alpinac_output_dir = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303_interesting_file\230311.0125.tank.1\alpinac_results\230311.0125.tank.1_EI_CI_output_198_mod")
    #path_for_alpinac_output_dir = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303_interesting_file\230311.0125.tank.1\alpinac_results\230311.0125.tank.1_EI_CI_output_198_mod")


    path_for_alpinac_input_dir = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303_interesting_file\230311.0125.tank.1\alpinac_results\230311.1025.tank.5\alpinac_input_files")
    path_for_alpinac_output_dir = Path(r"c:\Users\kaho\Desktop\data\data_Empa\results_campaign2023\Solvay_single_2\230311.1025.tank.5\alpinac_results\230311.1025.tank.5_EI_CI_output_4_mod")

    cmp_id = 4
    extended_quick_res_path = Path(r"C:\Users\kaho\Desktop\data\data_Empa\results_campaign2023\Solvay_single_2\230311.1025.tank.5")
    df_ext = pd.read_excel(r"c:\Users\kaho\Desktop\data\data_Empa\results_campaign2023\Solvay_single_2\230311.1025.tank.5\230311.1025.tank.5_peak_identification_results_extended.xlsx")
    create_modified_input_file(path_for_alpinac_input_dir, extended_quick_res_path, df_ext, cmp_id)

    get_alpinac_results(path_for_alpinac_output_dir, ionization_mode = "EI")


    db_dir = Path(r"C:\Users\kaho\Desktop\data\data_Empa\databases_Empa\TargetList")
    db_dir_extra_CH = Path(r"C:\Users\kaho\Desktop\data\data_Empa\databases_Empa\TargetList\nist_autodownloads\renamed_files")
    folder_path_campaign = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303_second_twin")


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

    result_dict = get_matchms_results_of_alpinac_fragments(path_for_alpinac_output_dir, dbs, df_peaks, rt_accept_int = (10, 10), ionization_mode = 'EI')