#%load_ext autoreload
#%autoreload
import glob
import itertools
import logging
import os
import chardet
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from alpinac_sups.read_h5 import H5_Reader
from alpinac.io_tools import MassCalibrationData
from alpinac.io_tools_hdf5 import hdf5Metadata
from alpinac.utils_data_extraction import f_2d_peak, extract_chomatograms
from db_reader import AlpinacFragmentsReader, JdxReader, AbstractReader
from matchms.Spectrum import Spectrum
from pathlib import Path
from scipy.signal import find_peaks
import pandas as pd
import matplotlib
from matchms.similarity import CosineGreedy
from matchms.filtering import normalize_intensities
from matchms.exporting import save_as_mgf
from pathlib import Path
from data_analysis_utils import formula_by_name, clean_empa_name_string
from alpinac.mass_calibration.utils import IonizationMethod



logging.getLogger("matchms").setLevel(logging.ERROR)


# defaul valued for the peak detection
# If rt in extracted file is smaller than peak list
TIME_FACTOR_NEG = 10
# if rt in extracted file is grater than in peak list
TIME_FACTOR_POS = 10
# Threshold to use for the RT PEAKS
RT_THRESHOLD = 1.5



# DICTIONARIES
# dictionary with manually added formulas
formulas_dict_manually_added = {
    'alpha-methylstyrene-b':['C2Cl3F3', 'prop-1-en-2-ylbenzene', '7407', 'InChI=1S/C9H10/c1-8(2)9-6-4-3-5-7-9/h3-7H,1-2H3', '2095.7'],
    'dichloro-trifluoro-pyridine':['C5Cl2F3N', '3,5-dichloro-2,4,6-trifluoropyridine', '61273', 'InChI=1S/C5Cl2F3N/c6-1-3(8)2(7)5(10)11-4(1)9', '2745.6'],
    #'CFC-113-only':['C2Cl3F3', '1,1,2-trichloro-1,2,2-trifluoroethane', '6428'],
    #'methoxyflurane_s':['C3H4Cl2F2O', '2,2-dichloro-1,1-difluoro-1-methoxyethane', '4119'],
    'C5H6O-a':['C5H6O', 'None', '10103117','None', '1856.7'],
    'C5H6O-b':['C5H6O', '2-methylfuran', '10797', 'InChI=1S/C5H6O/c1-5-3-2-4-6-5/h2-4H,1H3', '1865.9'],
    #'c-pentane':['C5H10', 'cyclopentane', '9253'],
    #'c-hexane':['C6H12', 'cyclohexane', '8078'],
    'DMS':['C2H6S', 'methylsulfanylmethane', '1068', 'InChI=1S/C2H6S/c1-3-2/h1-2H3', '1572'],
    #'1,2-dichloro-1,2-difluoroetheneE':['C2Cl2F2', '(E)-1,2-dichloro-1,2-difluoroethene', '3032334'], #trans/E isomerisch
    #'H-1211':['CBrClF2', 'bromo-chloro-difluoroethane', '9625'],
    #'1-chloro-1-propeneZ_s' :['C3H5Cl', '(Z)-1-chloro-1-propene', '5326315'],
    'ethylbenz':['C8H10', 'ethylbenzene', '7500', 'InChI=1S/C8H10/c1-2-8-6-4-3-5-7-8/h3-7H,2H2,1H3','2776.7'],
    'perfluorodimethyl-cyclohexane':['C8F16', '1,1,2,2,3,3,4,5,5,6-decafluoro-4,6-bis(trifluoromethyl)cyclohexane', '78975', 'InChI=1S/C8F16/c9-1(7(19,20)21)3(11,12)2(10,8(22,23)24)5(15,16)6(17,18)4(1,13)14]', '1779.1'],
    'perfluorodecaline':['C10F18', '1,1,2,2,3,3,4,4,4a,5,5,6,6,7,7,8,8,8a-octadecafluoronaphthalene', '9386', 'InChI=1S/C10F18/c11-1-2(12,5(17,18)9(25,26)7(21,22)3(1,13)14)6(19,20)10(27,28)8(23,24)4(1,15)16]', '2073.2'],
}


def prepare_extended_data_table(excel_quick_eval_file, formulas_dict_manually_added = {}, target_list = [], redo_analysis = False, formula_dict:dict = {}):
    """[summary]
    This function adds additional information to the excel file created by the quick_eval script (formulas, inchis, eval status, ...)
    """
    # actions:
    # sort by RT
    # add RT in min
    # add iupac name of cmp 
    # add formula of cmp
    # origin of cmp (NIST, mygu, xtra)
    # add difference from RT if available
    # add alpinac input file name for EI_CI and without N
    # redeinfe status in 10 categories dependent on score (0-0.1; 0.1-0.2, ..., 0.9-1)
    # add difference from RT if available
    file = excel_quick_eval_file
    # file from quick_eval (detect_peaks_script) #file = r"c:\Users\kaho\Desktop\data\data_Empa\Campaign202303_tests\230401.1747.air.3\230401.1747.air.3_peak_identification_results.xlsx"
    df = pd.read_excel(Path(file), sheet_name="Sheet1", index_col=0)
    print("excel file opened")
    # add info of results of extraction with mygu script
    # check if copy of excel file exists
    if os.path.isfile(file[:-5] + "_mod_kaho.xlsx"):
        print("modifed excel file already exists")
        if redo_analysis == True:
            # rename existing file to backup file
            os.rename(file[:-5] + "_mod_kaho.xlsx", file[:-5] + "_mod_kaho.xlsx_backup")
            print("backup file created")
        else: 
            df_copy = pd.read_excel(Path(file[:-5] + "_mod_kaho.xlsx"), sheet_name="Sheet1")
            print("modifed excel file opened")
    else:
        #print("modifed excel file created")
        # make copy of excel file
        df_copy = df.copy()
        # add index column to df_copy (ind_by_prom) at first position
        # get index of ascending intensity given in column "total_prominences" :
        ind_by_prom = np.argsort(np.array(-df["total_prominences"]))
        df_copy.insert(0, "ind_by_prom", ind_by_prom)
        # add column with RT in min after "RT"
        df_copy.insert(df.columns.get_loc("RT")+1, "RT_min", df_copy["RT"]/60)
        # add additional information
        for ind, matching_cmp in enumerate(["0_matching_cmp", "1_matching_cmp", "2_matching_cmp"]):
            identifier = [clean_empa_name_string(name_db) for name_db in df_copy[matching_cmp]]
            # generate list for results with same length as identifier
            formula = [""]*len(identifier)
            iupac_name = [""]*len(identifier)
            cid = [""]*len(identifier)
            inchi = [""]*len(identifier)
            RT_expected = [""]*len(identifier)
            # get formula, iupac_name, cid, inchi, RT_expected for each cmp
            for i, name_db in enumerate(identifier):

                formula_i, iupac_name_i, cid_i, inchi_i, RT_expected_i, formula_dict = formula_by_name(name = name_db, formula_dict = formula_dict, formulas_dict_manually_added  = formulas_dict_manually_added, target_list=target_list)
                if formula_i == "":
                    print(name_db)
                    print(formula_i, iupac_name_i, cid_i, inchi_i, RT_expected_i)
                formula[i] = formula_i
                iupac_name[i] = iupac_name_i
                cid[i] = cid_i
                inchi[i] = inchi_i
                RT_expected[i] = RT_expected_i

            #identif = [formula_by_name(name = clean_empa_name_string(name_db), formula_dict = formula_dict, formulas_dict_manually_added  = formulas_dict_manually_added, target_list=target_list) for name_db in df_copy[matching_cmp]]
            #formula, iupac_name, cid, inchi, RT_expected, formula_dict = zip(*identif)   
            # get column index of 0_matching_cmp
            col_ind_cmp = df_copy.columns.get_loc(matching_cmp)
            # insert the lines after the 0_matching_cmp column
            df_copy.insert(col_ind_cmp+1, str(ind) + "_name", iupac_name)
            df_copy.insert(col_ind_cmp+2, str(ind) + "_formula", formula)
            df_copy.insert(col_ind_cmp+3, str(ind) + "_inchi", inchi)
            df_copy.insert(col_ind_cmp+4, str(ind) + "_RT_expected", RT_expected)

        # print out the following columns  0_matching_cmp, 0_RT_expected
        print(df_copy[["0_matching_cmp", "0_RT_expected"]])
        # add column with alpinac input file name for EI_CI and without N
        df_copy["alpinac_input_file_EI_CI"] = [str(Path(file).stem.split("_")[0] + "_EI_CI_input_")+str(cmp_id)+".txt" for cmp_id in df_copy["cmp_id"]]
        #df_copy["alpinac_input_file_EI_CI_no_N"] = [str(Path(file).stem.split("_")[0] + "EI_CI_no_N_input_")+str(cmp_id)+".txt" for cmp_id in df_copy["cmp_id"]]

        # add empty column
        df_copy["extraction_status"] = 0 # Extraction: 0 = not run, 1 = extracted, 2 = extracted, but failed, 3 = extracted, but too litte peaks
        df_copy["extraction_number_cmps"] = "" 
        df_copy["extraction_peaks_per_cmp_EI"] = ""
        df_copy["extraction_peaks_per_cmp_CI"] = ""
        df_copy["alpinac_status"] = 0 # Alpinac: 0 = not run, 1 = run, 2 = run, but failed
        df_copy["alpinac_EI"] = ""
        df_copy["alpinac_CI"] = ""
        df_copy["alpinac_no_N_status"] = 0 # Alpinac: 0 = not run, 1 = run, 2 = run, but failed
        df_copy["alpinac_EI_no_N"] = ""
        df_copy["alpinac_CI_no_N"] = ""
        # save as new excel file
        #df_copy.to_excel(file[:-5] + "_mod_kaho.xlsx", index=False)
        #sort df by retention time
        df_copy2 = df_copy.sort_values(by=["RT"])

        # save df_copy as excel file
        df_copy2.to_excel(file[:-5] + "_extended.xlsx", index=False)
    return df_copy2, formula_dict

def creates_results_df(
        compounds: list[int], 
        extracted_rts_peak: np.ndarray, 
        extracted_rts_start: np.ndarray, 
        extracted_rts_end: np.ndarray,
        number_of_peaks: np.ndarray,
        total_intesities: np.ndarray,
        total_prominences: np.ndarray,
        mask_identified: np.ndarray,
        mask_possibly_identified: np.ndarray,
        scores: dict[int, dict[str, float]],
        n_firsts: int=150) -> pd.DataFrame:
    """
    Parameters
    ----------
    compounds : list[int]
        List of compounds to consider
    extracted_rts_peak : np.ndarray
        RT of the peak
    extracted_rts_start : np.ndarray
        RT of the start of the peak
    extracted_rts_end : np.ndarray
        RT of the end of the peak
    number_of_peaks : np.ndarray
        Number of peaks in the compound
    total_intesities : np.ndarray
        Total intensity of the peaks
    total_prominences : np.ndarray
        Total prominence of the peaks
    mask_identified : np.ndarray
        Mask of the identified compounds
    mask_possibly_identified : np.ndarray
        Mask of the possibly identified compounds
    scores : dict[int, dict[str, float]]
        Scores of the compounds
    n_firsts : int, optional
        Number of compounds to consider, by default n_firsts
    """
    missing_compounds = []
    df = pd.DataFrame(index=range(n_firsts), columns=[
        "cmp_id",
        "RT",
        "rt_start",
        "rt_end",
        "number_of_peaks",
        "total_intesities",
        "total_prominences",
        "Status",
        "0_matching_cmp",
        "0_score",
        "1_matching_cmp",
        "1_score",
        "2_matching_cmp",
        "2_score",
    ])
    for ranking, cmp_id in enumerate(compounds[:n_firsts]):
        line = {
            'cmp_id': cmp_id,
            'RT': extracted_rts_peak[cmp_id],
            'rt_start': extracted_rts_start[cmp_id],
            'rt_end': extracted_rts_end[cmp_id],
            'number_of_peaks': number_of_peaks[cmp_id],
            'total_intesities': total_intesities[cmp_id],
            'total_prominences': total_prominences[cmp_id],
        }
        #print(
        #    cmp_id,
        #    f"peak:{extracted_rts_peak[cmp_id]:.2f}±{extracted_rts_peak_std[cmp_id]:.2f},start:{extracted_rts_start[cmp_id]:.2f}, end:{extracted_rts_end[cmp_id]:.2f}",
        #)
        if mask_identified[cmp_id]:
            #print("Was  identified", matched_cmps[cmp_id])
            #if almost_matched_cmps[cmp_id]:
                #print("Could also have been ", almost_matched_cmps[cmp_id])
            line['Status'] = 'Identified'
        else:
            if mask_possibly_identified[cmp_id]:
                #print("Was almost identified", almost_matched_cmps[cmp_id])
                line['Status'] = 'Almost Identified'
            else:
                #print("Was not identified")
                line['Status'] = 'Unknown'
            missing_compounds.append(cmp_id)
        # Get the score for that compound by score
        cmp_names = list(scores[cmp_id].keys())
        cmp_scores = list(scores[cmp_id].values())

        for i, index in enumerate(np.argsort(cmp_scores)[::-1]):
            if i == 3:
                break
            line[f"{i}_matching_cmp"] = cmp_names[index]
            line[f"{i}_score"] = cmp_scores[index]

        df.loc[ranking] = line

    return df

def detect_peaks(
        file: Path,
        rt_start: float,
        rt_end: float,
        ionisaztion_dict: dict,
        db_dir: Path,
        db_dir_extra: Path,
        use_substances_with_unknown_RT: bool,
        out_dir: Path = None,
        inoization_to_use: str = 'EI',
        d_m: float = 0.1,
        masses_to_ignore: list = [185, 187],
        n_firsts: int = 150,
        time_factor_neg: int = TIME_FACTOR_NEG,
        time_factor_pos: int = TIME_FACTOR_POS,
        rt_threshold: float = RT_THRESHOLD):
    """
    Parameters
    ----------
    file : Path
        Path to the h5 file
    rt_start : float
        Start of the RT to consider
    rt_end : float
        End of the RT to consider
    ionisaztion_dict : dict
        Dictionary of the ionization segments to use for different ionization modes
    db_dir : Path
        Path to the database
    db_dir_extra_CH : Path
        Path to the extra database for extra libraries. These were not in the Empa confirmed db (with RT) but additional ones saved with default NIST names (CAS)
    use_substances_with_unknown_RT : bool
        If True, also use substances with unknown RT
    out_dir : Path
        Path to the output directory
    inoization_to_use : str
        Ionization mode to use
    d_m : float
        Mass window to use for the peak detection
    masses_to_ignore : list
        List of masses to ignore
    n_firsts : int
        Number of RT peaks to consider
    TIME_FACTOR_NEG : int
        Factor to use for the negative time factor
    TIME_FACTOR_POS : int
        Factor to use for the positive time factor
    RT_THRESHOLD : float
        Threshold to use for the RT
    """

    logger = logging.getLogger("detect_peaks")
    logging.basicConfig(level=logging.INFO)     


    if out_dir is None:
        out_dir = file.with_suffix('')
    out_dir.mkdir(exist_ok=True)

    # get all hdf5 files in the directory

    # Read the data from the h5 file
    # h5_reader = H5_Reader(file)
    h5_reader = H5_Reader(file, ionisaztion_dict)
    h5_helper = hdf5Metadata(file, mode="r")

    # find for each substances in the known peaks the values of their masses

    df_peaks = pd.read_csv(db_dir / "TargetList.csv")
    # kaho 230622: add to do a pseudo target list for CH combis, with no known RT, evenutally manually limit RTs later
    update_pseudo_peaks = False
    if update_pseudo_peaks:
        xtra_CH_db = JdxReader(db_dir_extra_CH).read_all()
        df_pseudo_peaks = []
        # try to read it if it exists
        try:
            df_pseudo_peaks = pd.read_csv(db_dir / "PseudoTargetList.csv", names = df_peaks.columns)
        except:
            pass
        # get the keys of xtra_CH_dbwhich are not in the df_pseudo_peaks
        if len(df_pseudo_peaks) == 0:
            # set names of columns
            df_pseudo_peaks = pd.DataFrame(columns = df_peaks.columns)
            df_pseudo_peaks['Substance'] = xtra_CH_db.keys()
            df_pseudo_peaks['RT'] = 1700
            df_pseudo_peaks['RT-Window'] = 1699
            for i in range(len(df_pseudo_peaks)):  
                # get largest peaks
                intensities = xtra_CH_db[df_pseudo_peaks['Substance'][0]].peaks.intensities
                # get masses belonging to the three largest peaks
                masses = xtra_CH_db[df_pseudo_peaks['Substance'][0]].peaks.mz
                df_pseudo_peaks['Target'], df_pseudo_peaks['Q1'],  df_pseudo_peaks['Q2']= masses[np.argsort(intensities)[::-1][:3]]  
            df_pseudo_peaks.to_csv(db_dir / "PseudoTargetList.csv", index = False)
        else: 
            raise NotImplementedError("PseudoTargetList.csv already exists, please check if it is correct and delete it if not")
        # Add pseudo peaks to the list
        # save the pseudo peaks
        
    df_pseudo_peaks = pd.read_csv(db_dir / "PseudoTargetList.csv")
    if use_substances_with_unknown_RT:
        df_peaks = pd.concat([df_peaks, df_pseudo_peaks], ignore_index=True)


    # Load the TOF data
    # full_tof_data = h5_reader.read_full_tof_data()
    # full_tof_data = full_tof_data.reshape(full_tof_data.shape[0], full_tof_data.shape[2])
    full_tof_data = h5_reader.read_full_ionization_data(inoization_to_use)


    #  Find where are the retention times needed and what is the data we want to extract
    rt_id_start, rt_id_end = np.searchsorted(h5_reader.time_axis, (rt_start, rt_end))
    tof_data = full_tof_data[rt_id_start:rt_id_end, :]
    time_axis = h5_reader.time_axis[rt_id_start:rt_id_end]

    # Small function to convert from bin to rt
    time_bin_to_time = lambda x: np.interp(x, np.arange(len(time_axis)), time_axis)

    masses = h5_reader.get_mass_axis_of(rt_start, rt_end, spectra_mode = inoization_to_use)

    #Extraction of the peaks from the data

    peaks_bins = {}
    all_peak_params = {}


    chromato_dict = extract_chomatograms(
        tof_data,
        masses,
        #compounds
        masses=[i for i in range(2, 300) if i not in masses_to_ignore],
        d_m=d_m
    )

    for mass, chromato in chromato_dict.items():
        # Find the peaks on each chomatorogram
        peaks_timebin, peaks_params = find_peaks(
            chromato,
            # Filter out data when peaks are too small
            width=3,
            height=max(10, np.median(chromato) ),
        )
        peaks_bins[mass] = peaks_timebin
        all_peak_params[mass] = peaks_params
    p, n_peaks = np.unique(
        np.concatenate([peaks for peaks in peaks_bins.values()]), return_counts=True
    )
    # Show the peaks
    plt.figure()
    for m, peaks in peaks_bins.items():
        # plt.scatter(m * np.ones_like(peaks), time_axis[peaks], s=3*np.log(all_peak_params[m]['prominences']))
        left = time_bin_to_time(all_peak_params[m]["left_ips"])
        right = time_bin_to_time(all_peak_params[m]["right_ips"])
        plt.bar(m * np.ones_like(peaks), bottom=left, height=right - left, width=2 * d_m)

    # prepare the data for the assignement of peaks

    # merge all peak data in single arrays 
    merged_peaks_params = {
        key: np.concatenate([all_peak_params[mass][key] for mass in peaks_bins.keys()])
        for key in peaks_params.keys()
    }
    # Convert bins to RTs
    # Left is where the peak starts and ends on the right
    merged_peaks_params['start_rt'] = time_bin_to_time(
        merged_peaks_params["left_ips"] 
    )
    merged_peaks_params['end_rt'] = time_bin_to_time(
        merged_peaks_params["right_ips"] 
    )
    merged_peaks_params['peak_rt']  = time_bin_to_time(
        np.concatenate([peaks_bins[mass] for mass in peaks_bins.keys()])
    )

    # get the mass for which the peak was found
    merged_peaks_params['mass']  = np.concatenate(
        [np.full(len(peaks), mass) for mass, peaks in peaks_bins.items()]
    )

    n_peaks = len(merged_peaks_params['start_rt'])

    # Assing peaks to compounds
    cmp_id = 0

    mask_assigned = np.zeros(n_peaks, dtype=bool)
    cmp_ids = np.empty(n_peaks, dtype=int)

    end_rt, start_rt = merged_peaks_params['end_rt'], merged_peaks_params['start_rt']
    dist_rt = end_rt - start_rt

    while mask_assigned.sum() < n_peaks:
        # Take the first not assigned peak
        # TODO: define what strategy would be better to start
        non_assigned_peaks = np.argwhere(~mask_assigned)[0]
        # Take the smallest peak so that noises assigned can be corrected at the end by bigger peaks
        arg_min = np.argmin(merged_peaks_params['prominences'][non_assigned_peaks])
        peak_index = non_assigned_peaks[arg_min]
        # print(f"Peak {peak_index} is assigned to {cmp_id=}")

        mask_this_compound = (
            (
                # Start of the peak is quite close to the RT
                (start_rt > start_rt[peak_index] - rt_threshold)
                & (start_rt < start_rt[peak_index] + rt_threshold)
            )
            & (
                # End of the peak is quite close to the end RT
                (end_rt > end_rt[peak_index] - rt_threshold)
                & (end_rt < end_rt[peak_index] + rt_threshold)
            )
            # & (
            #     # Distance (peak width) is similar
            #     (dist_rt > dist_rt[peak_index] - RT_THRESHOLD)
            #     & (dist_rt < dist_rt[peak_index] + RT_THRESHOLD)
            # )
        )
        # peaks in this compound
        cmp_ids[mask_this_compound] = cmp_id
        cmp_id += 1
        mask_assigned[mask_this_compound] = True

    n_compounds = cmp_id


    # Some times the peaks will be very close to anothers
    # TODO: check that the masses present are different
    # TODO: check that the rt are close to each other

    # Load the data bases
    dbs: dict[str, dict[str, Spectrum]] = {
        "nist": JdxReader(db_dir / "nist_db").read_all(),
        "myriam": AlpinacFragmentsReader(db_dir / "myriam_ei").read_all(),
        "xtra_CH": JdxReader(db_dir_extra_CH).read_all(),
    }

    # Analyse the compounds now and try to make then match the database

    cosine_greedy = CosineGreedy(tolerance=0.2)
    matches = {}
    scores: dict[int, dict[str, float]] = {}
    mask_identified = np.zeros(n_compounds, dtype=bool)
    mask_possibly_identified = np.zeros(n_compounds, dtype=bool)
    total_intesities = np.zeros(n_compounds, dtype=float)
    total_prominences = np.zeros(n_compounds, dtype=float)
    extracted_rts_start = np.zeros(n_compounds, dtype=float)
    extracted_rts_end = np.zeros(n_compounds, dtype=float)
    extracted_rts_peak = np.full(n_compounds, np.nan, dtype=float)
    extracted_rts_peak_std = np.full(n_compounds, np.nan, dtype=float)
    number_of_peaks = np.full(n_compounds, np.nan, dtype=float)
    matched_cmps = [[] for _ in range(n_compounds)]
    almost_matched_cmps = [[] for _ in range(n_compounds)]
    ms_spectras = []
    spectra_dir = out_dir / "matchms_spectra"
    spectra_dir.mkdir(exist_ok=True)


    for cmp_id in range(n_compounds):
        mask_cmp = cmp_ids == cmp_id
        n_peaks = (mask_cmp).sum()
        if n_peaks < 2:
            continue

        cmp_start_rt = start_rt[mask_cmp].mean()
        cmp_end_rt = end_rt[mask_cmp].mean()
        extracted_rts_start[cmp_id] = cmp_start_rt
        extracted_rts_end[cmp_id] = cmp_end_rt

        extracted_rts_peak[cmp_id] = merged_peaks_params["peak_rt"][mask_cmp].mean()
        extracted_rts_peak_std[cmp_id] = merged_peaks_params["peak_rt"][mask_cmp].std()

        mz = merged_peaks_params["mass"][mask_cmp]
        i = merged_peaks_params["peak_heights"][mask_cmp]
        number_of_peaks[cmp_id] = np.sum(mask_cmp)
        total_intesities[cmp_id] = np.sum(i)
        total_prominences[cmp_id] = np.sum(merged_peaks_params["prominences"][mask_cmp])

        # Get compounds matching with the rt
        mask_matching_time = (
            cmp_start_rt > df_peaks["RT"] - time_factor_neg * df_peaks["RT-Window"] / 2
        ) & (cmp_end_rt < df_peaks["RT"] + time_factor_pos * df_peaks["RT-Window"] / 2)
        compounds = list(df_peaks.loc[mask_matching_time, "Substance"])

        matchms_spectra = Spectrum(mz=np.array(mz, dtype=float), intensities=i)


        with open(spectra_dir / f"{cmp_id}.mgf", 'w') as f:
            save_as_mgf([matchms_spectra], f)
        ms_spectras.append(matchms_spectra)
        # Compare the compound to known compounds

        scores[cmp_id] = {}
        for compound in compounds:

            for db, db_compounds in dbs.items():

                if compound not in db_compounds:
                    #print("Unknown RT of", compound, "in db", db)
                    continue
                # Compare the two spectra
                reference = db_compounds[compound]
                # use only masses in the range of the experimental device
                mask_reference_in_range = (reference.peaks.mz > 27) & (reference.peaks.mz < 500)
                reference_in_range  = Spectrum(mz=reference.peaks.mz[mask_reference_in_range], intensities=reference.peaks.intensities[mask_reference_in_range], metadata={"precursor_mz": -1})
                score_res = cosine_greedy.pair(
                    normalize_intensities(reference_in_range), normalize_intensities(matchms_spectra)
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
                        logger.warning(f"More peaks matched than expected in measured spectra (mz>27, peak hight 5% of max), adapt measured range, in the reference for {compound} in {db}")
                    else:
                        score = w_cossim*score_cossim + w_no_peaks * matched_peaks/total_no_peaks_in_lim
                # Add score to the results dictionaries
                scores[cmp_id][f"{compound}_{db}"] = score
                matches[compound][cmp_id] = score
                # Decide what to do based on the score
                if score > 0.9:
                    if mask_identified[cmp_id]:
                        print(cmp_id, "Was already identified")
                    mask_identified[cmp_id] = True
                    matched_cmps[cmp_id].append(compound)
                elif score > 0.1:
                    mask_possibly_identified[cmp_id] = True
                    almost_matched_cmps[cmp_id].append(compound)
                else:
                    pass
        # Set attributed or not based on the score

    # find the largest compuond remaining
    cmp_by_total_i = np.argsort(total_prominences)[::-1]
    df = creates_results_df(compounds=cmp_by_total_i, 
                            extracted_rts_peak=extracted_rts_peak, 
                            extracted_rts_start=extracted_rts_start, 
                            extracted_rts_end=extracted_rts_end, 
                            number_of_peaks=number_of_peaks, 
                            total_intesities=total_intesities, 
                            total_prominences=total_prominences, 
                            mask_identified=mask_identified, 
                            mask_possibly_identified=mask_possibly_identified, 
                            scores=scores, 
                            n_firsts=n_firsts)
    df.to_excel(
        out_dir / (file.stem + "_peak_identification_results.xlsx"),
    )
    df

    #Saves a peak list from the identified peaks 
    df_peaklist = pd.DataFrame({
        "RT": extracted_rts_peak[mask_identified],
        "compound": [matched_cmps[i][0] if len(np.unique(matched_cmps[i])) == 1 else matched_cmps[i] for i in np.argwhere(mask_identified).reshape(-1)],}
    )
    df_peaklist

    return df_peaklist, df, df_peaks, matched_cmps, almost_matched_cmps, ms_spectras, extracted_rts_start, extracted_rts_end, extracted_rts_peak


# main
if __name__ == "__main__":


    print(formula_by_name('463-49-0',{}))
    print(formula_by_name('3-ethylpentane',{}))
    print(formula_by_name('617-78-7',{}))



    # Where the library can be found
    db_dir = Path(r"C:\Users\kaho\Desktop\data\data_Empa\databases_Empa\TargetList")
    db_dir_extra_CH = Path(r"C:\Users\kaho\Desktop\data\data_Empa\databases_Empa\TargetList\nist_autodownloads\renamed_files")
    folder_path_campaign = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303")
    folder_path_campaign = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303_interesting_file")
    folder_path_campaign = Path(r"G:\503_Themen\Klima\kaho\results_campaign2023\standard2")

    path_target_list = db_dir/Path("TargetList_extended.csv")

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



    inoization_to_use = "EI"
    ionisaztion_dict = {
        # Last value is not included
        "CI": (1, 2),  # Segment 1
        # "EI": (0, 1),  # Segment 0 (For old files with only EI)
        "EI": (3, 4),
    }
    
    rt_start, rt_end = 1400, 3050
    d_m = 0.1
    # These are some masses that we don't want as they have maybe a column bleed
    masses_to_ignore = [185, 187]

    h5_files = [f for f in folder_path_campaign.glob("**/*.h5") if f.with_suffix('.EI_mc.txt').exists()]


    logger = logging.getLogger("detect_peaks")
    logging.basicConfig(level=logging.INFO)

    for file in h5_files:
        print(file)
        filename_quick_results = file.with_suffix('')/(Path(file.name).stem + "_peak_identification_results.xlsx")
        # check if alrady processed
        if filename_quick_results.exists():
            print("already processed")
            continue

        df_peaklist, df, df_peaks, matched_cmps, almost_matched_cmps, ms_spectras, extracted_rts_start, extracted_rts_end, extracted_rts_peak = detect_peaks(
            file = file,
            rt_start = rt_start,
            rt_end = rt_end,
            ionisaztion_dict = ionisaztion_dict,
            db_dir = db_dir,
            db_dir_extra = db_dir_extra_CH,
            use_substances_with_unknown_RT=True,
            out_dir=file.with_suffix(''),
            inoization_to_use=inoization_to_use,
            d_m=0.1,
            masses_to_ignore=[185, 187],
            n_firsts=150,
            time_factor_neg=TIME_FACTOR_NEG,
            time_factor_pos=TIME_FACTOR_POS,
            rt_threshold=RT_THRESHOLD,
        )


    optional_specific_plot = False
    cmp_id = 1124

    if optional_specific_plot == True:
        # Optianally show some specific peaks
        # access the mising peaks data

        extracted_rts_start[cmp_id], extracted_rts_end[cmp_id], extracted_rts_peak[cmp_id], almost_matched_cmps[cmp_id], ms_spectras[cmp_id], 
        ms_spectras[cmp_id].plot()


    # GET UPDATED TARGET LIST
    # figure out encoding of target list
    with open(path_target_list, 'rb') as f:  
        result = chardet.detect(f.read())  # or readline if the file is large
    target_list_encoding = result['encoding']

    # TARGET LIST
    # load the target list
    target_list = pd.read_csv(path_target_list, encoding=target_list_encoding)
    update_target_list = False
    # update_target_list if necessary:
    if update_target_list:
        target_list_iupac = []
        target_list_inchi = []
        # add a column with iupac_name and inchi_key
        for subs in target_list.Substance:
            iupac_name = ''
            inchi = ''
            try:
                formula, iupac_name , _, inchi,_, formulas_dict = formula_by_name(name=subs, formulas_dict_manually_added=formulas_dict_manually_added, formula_dict = formula_dict, target_list=target_list)
            except:
                pass
            target_list_iupac.append(iupac_name)
            target_list_inchi.append(inchi)
        target_list['iupac_name'] = target_list_iupac
        target_list['inchi_key'] = target_list_inchi
        target_list.to_csv(path_target_list, index=False)
        print('Target list updated')

    # CREATE EXTENDED DATA TABLES FOR ALL EXCEL FILES

    dir_contents = os.listdir(folder_path_campaign)
    all_excel_files = []
    all_excel_files_extended = []
    extension1 = "*results.xlsx"
    extension2 = "*results_extended.xlsx"
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
                # print the list of excel files
                list_excel.extend(excel_files)

    # do prepare_extended_data_table for all excel files
    for file in all_excel_files:
        # check if file already exists  
        if os.path.isfile(file[:-5] + "_extended.xlsx"):
            print("extended excel file already exists")
        else:
            print("processing file: " + file)
            df_ext, formula_dict = prepare_extended_data_table(excel_quick_eval_file=file, formulas_dict_manually_added=formulas_dict_manually_added, target_list=target_list, redo_analysis=False, formula_dict=formula_dict)
            # extend all_excel_files_extended
            all_excel_files_extended.append(file[:-5] + "_extended.xlsx")

    # update formulas_dict
    with open( db_dir.parent /Path("target_list_supplementary_files")/ 'formula_dict.txt', 'w') as f:
            for key, val in formula_dict.items():
                f.write('%s:%s %s' % (key, val, '\n')) # python will convert \n to os.linesep   

