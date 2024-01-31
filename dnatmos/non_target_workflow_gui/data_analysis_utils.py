import os
import time
from pathlib import Path
from typing import Union
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

from matplotlib import pyplot as plt
from matplotlib.lines import segment_hits
import numpy as np
from mendeleev import element
from scipy.signal import find_peaks, find_peaks_cwt


from math import isclose
from collections import Counter
from typing import List
# import chem_formula
from molmass import Formula as chem_formula
from alpinac_sups.read_h5 import read_desc, read_ionisation_modes, read_mass_axis, H5_Reader, read_time_axis, read_mass_cal_par
from scipy.optimize import minimize
#import StickSpectra class
from alpinac_gui.matchms_funcs import StickSpectra
from pyvalem.formula import Formula
import pubchempy as pcp
import cirpy
from alpinac_gui.matchms_funcs import formula_from_namestring, formula_to_latex, halocarbname2formula, plot_mass_spectra, load_from_mgf, metadata_processing, peak_processing, AlpinacData, halocarbname2formula



# FUNCTIONS 
# def index_to_mass(i, p1, p2):
#     return ((np.array(i)-p2)/p1)**2

# def mass_to_index(mass, p1, p2):
#     return np.sqrt(mass*p1 + p2)

# def rt_index_to_rt(i, s_per_bin, offset):
#     return (np.array(i) + offset)*s_per_bin

# def rt_to_rt_index(rt, s_per_bin, offset):
#     return (np.array(rt)/s_per_bin) - offset

def mass_index_to_mass(mass_index, params:list = [2835.23579263574, -4049.42781522328], mode: int = 0) -> np.ndarray:
    """ This function converts the mass_index to the mass value.
    Parameters: 
    mass_index: index of the mass value
    params: list of parameters for the conversion. The first parameter is the number of seconds per bin, the second parameter is the offset.
    mode: mode of the conversion. 0: square root conversion"""
    if mode == 0:
        p1 = params[0]
        p2 = params[1]
        return_array = ((np.array(mass_index)-p2)/p1)**2
        return_func = lambda x: ((x-p2)/p1)**2
    elif mode == 2:
        p1 = params[0]
        p2 = params[1]
        p3 = params[2]
        return_array = ((np.array(mass_index)-p2)/p1)**(1/p3)
        return_func = lambda x: ((x-p2)/p1)**(1/p3)
    else:
        print("Error: Invalid mode. Choose from 0 or 2")
        return None, None
    return return_array, return_func

def mass_to_mass_index(mass, params: list = [2835.23579263574, -4049.42781522328], mode: int = 0) -> np.ndarray:
    """ This function converts the mass value to the mass_index.
    Parameters:
    mass: mass value
    params: list of parameters for the conversion. The first parameter is the number of seconds per bin, the second parameter is the offset.
    mode: mode of the conversion. 0: square root conversion
    Returns:
    mass_index: index of the mass value
    """
    if mode == 0:
        p1 = params[0]
        p2 = params[1]
        return_array = np.sqrt(mass*p1 + p2)
        return_func = lambda x: np.sqrt(x*p1 + p2)
    elif mode == 2:
        p1 = params[0]
        p2 = params[1]
        p3 = params[2]
        return_array = (mass*p1 + p2)**p3
        return_func = lambda x: (x*p1 + p2)**p3
    else:
        print("Error: Invalid mode. Choose from 0 or 2")
        return None, None
    return return_array, return_func    

def rt_index_to_rt(rt_index, params:list = [0.2,0], mode:int = 0) -> np.ndarray:
    """ This function converts the rt_index to the rt value.
    Parameters:
    rt_index: index of the rt value
    params: list of parameters for the conversion. The first parameter is the number of seconds per bin, the second parameter is the offset.
    mode: mode of the conversion. 0: linear conversion
    Returns:
    rt: rt value
    """
    if mode == 0:
        s_per_bin = params[0]
        offset = params[1]
        return_array = np.array(rt_index)*s_per_bin + offset
        return_func = lambda x: x*s_per_bin + offset
    return return_array, return_func    

def rt_to_rt_index(rt, params:list = [0.2,0], mode:int = 0) -> np.ndarray:
    """ This function converts the rt value to the rt_index.
    Parameters:
    rt: rt value
    params: list of parameters for the conversion. The first parameter is the number of seconds per bin, the second parameter is the offset.
    mode: mode of the conversion. 0: linear conversion
    Returns:
    rt_index: index of the rt value
    """
    if mode == 0:
        s_per_bin = params[0]
        offset = params[1]
        return_array = (np.array(rt)/s_per_bin) - offset
        return_func = lambda x: (x/s_per_bin) - offset
    return return_array, return_func

def load_data(file_path: Union[str, Path], ionisation_dict, ionisation_type: str) -> Union[None, tuple]:
    """
    This function loads the data from the given file_path and returns the data and the reader object.
    Parameters:
    file_path: path to the file
    ionisation_dict: dictionary with the ionisation types as keys and the corresponding segment numbers as values
    ionisation_type: ionisation type to be loaded. Choose from 'EI' or 'CI'
    Returns:
    data: data as numpy array
    reader: reader object  
    """
    try:
        file = Path(file_path)
        
        if not file.exists():
            raise FileNotFoundError("File not found")
        
        if ionisation_type not in ionisation_dict:
            raise ValueError("Invalid ionisation type. Choose from 'EI' or 'CI'")
        
        reader = H5_Reader(file, ionisation_dict)

        
        if ionisation_type == 'EI':
            data = reader.read_full_ionization_data('EI')
        elif ionisation_type == 'CI':
            data = reader.read_full_ionization_data('CI')
        
        return data, reader
    
    except FileNotFoundError as e:
        print(f"Error: {e}")
        return None
    
    except ValueError as e:
        print(f"Error: {e}")
        return None
    
    except Exception as e:
        print(f"Error: {e}")
        return None
    
def get_sum(data, sum_type, range_starts:np.array = None, range_ends: np.array =None, tol=0.1):
    """ This function returns the sum of the data in the given range(s) of the given sum_type
    sum_type: "RT" (for RT_index spectrum) or "mass" (for mass index spectrum)
    range_starts: list of start values of the ranges (in index units)
    range_ends: list of end values of the ranges (index untis)
    tol: tolerance for the range values
    """
    " range_starts and range_ends must be lists of the same length"
    #if range_start is not a list, make it a list
    # test if range_starts is a numpy array
    range_starts = np.atleast_1d(range_starts)
    range_ends = np.atleast_1d(range_ends)

    if range_starts is not None and not isinstance(range_starts, np.ndarray):
        range_starts = np.array(range_starts)
    #if range_end is not a list, make it a list
    if range_ends is not None and not isinstance(range_ends, np.ndarray):
        range_ends = np.array(range_ends)
    if sum_type == "RT":
        if range_starts[0] is not None and range_ends[0] is not None:
            data_sum = np.zeros((data.shape[0], len(range_starts)))
            for i, (range_start, range_end) in enumerate(zip(range_starts, range_ends)):
                rt_values = np.arange(data.shape[1])
                range_mask = np.logical_and(rt_values >= range_start-tol, rt_values <= range_end+tol)
                data_sum[:, i] = data[:, range_mask].sum(axis=1)
        else:
            data_sum = data.sum(axis=1)
    elif sum_type == "mass": # 
        if range_starts[0] is not None and range_ends[0] is not None:
            data_sum = np.zeros((data.shape[1], len(range_starts)))
            for i, (range_start, range_end) in enumerate(zip(range_starts, range_ends)):
                print (range_start, range_end)
                # plot data for i
                #fig, ax = plt.subplots(1,1)
                # assign index to mass
                #index = np.arange(data.shape[1])
                #ax.plot(index_to_mass(index, p1=p1, p2=p2), data[650:750,:].sum(axis=0))
                mass_values = np.arange(data.shape[0])
                range_mask = np.logical_and(mass_values >= range_start-tol, mass_values <= range_end+tol)
                data_sum[:, i] = data[range_mask, :].sum(axis=0)
                #ax.plot(index_to_mass(index, p1=p1, p2=p2), data_sum[i, :])
        else:
            # if no ranges are given, sum over all RT values
            data_sum = data.sum(axis=0)
    else:
        print("Error: Invalid sum type. Choose from 'RT' or 'mass'")
        return None
    # if only one range is given, return the sum as a 1D array
    if data_sum.shape[0] == 1:
        data_sum = data_sum[0, :]
    return data_sum

def get_mass_spec_for_RT_range(data, ranges_start:np.array = None, ranges_end:np.array = None, rt_index_to_rt_mode:int = 0, rt_index_to_rt_params:list = [float(1/5),0], mass_index_to_mass_mode: int = 0, mass_index_to_mass_params:list = [2835.23579263574, -4049.42781522328], range_tol: float=0.1, range_mode:str = 'abs', sum_x:str = 'mass', normalize:bool = False):
    """
    This function returns the mass spec for the given RT range{s}.
    Parameters:
    data: GC-MS data
    ranges_start: list of start values of the ranges (in RT or index units, depending on the range_mode)
    ranges_end: list of end values of the ranges (RT units or index units, depending on the range_mode)
    rt_index_to_rt_mode: mode of the conversion from RT index to RT value. 0: linear conversion
    rt_index_to_rt_params: list of parameters for the conversion from RT index to RT value. The first parameter is the number of seconds per bin, the second parameter is the offset.
    mass_index_to_mass_mode: mode of the conversion from mass index to mass value. 0: square root conversion
    mass_index_to_mass_params: list of parameters for the conversion from mass index to mass value. The first parameter is the number of seconds per bin, the second parameter is the offset.
    range_tol: tolerance for the range values (in RT or index units, depending on the range_mode)
    range_mode: mode of the range values. 'abs': absolute values, 'ind': index values
    sum_x: x axis of the sum. 'mass': mass values, 'mass_index': mass index values
    normalize: if True, normalize the mass spec
    Returns:
    mass_spec_df: pandas dataframe with the mass spec data. The columns are 'mass' and 'intensity' (or 'intensity_norm' if normalize is True) 
    """
    # get the mass spec for the given RT range
    if range_mode == 'abs' and ranges_start is not None and ranges_end is not None: # range is given in absolute values not in indici
        # convert RTs to RT indici
        ranges_start, _ = rt_to_rt_index(ranges_start, params=rt_index_to_rt_params, mode=rt_index_to_rt_mode)
        ranges_end, _ = rt_to_rt_index(ranges_end, params=rt_index_to_rt_params, mode=rt_index_to_rt_mode)
        range_tol, _ = rt_to_rt_index(range_tol, params=rt_index_to_rt_params, mode=rt_index_to_rt_mode)
        # convert to int
        ranges_start = ranges_start.astype(int)
        ranges_end = ranges_end.astype(int)
        #range_tol = range_tol.astype(int)
    # get the sum of the data in the given range
    mass_spec = get_sum(data, sum_type="mass", range_starts=ranges_start, range_ends=ranges_end, tol=range_tol)
    # convert the mass indici to mass values
    mass_index = np.arange(mass_spec.shape[0])  
    if sum_x == 'mass':
        mass_x, _ = mass_index_to_mass(mass_index, params=mass_index_to_mass_params, mode=mass_index_to_mass_mode)
    elif sum_x == 'mass_index':
        mass_x = mass_index
    else:
        print("Error: Invalid sum_x. Choose from 'mass' or 'mass_index'")
        return None

    mass_spec_df = pd.DataFrame()
    mass_spec_df['mass'] = mass_x
    mass_spec_df['intensity'] = mass_spec

    # convert to pandas dataframe
    # generate data
    # data_mass_spec = np.array([mass_x, mass_spec]).T    
    # mass_spec_df = pd.DataFrame(data_mass_spec, columns=[sum_x, 'intensity'])
    if normalize:
        #mass_spec_df = mass_spec_df.rename(columns={'intensity': 'intensity_norm'})
        mass_spec_df['intensity_norm'] = mass_spec/mass_spec.sum()
    # convert to ContSpec object # TODO: write such a class!
    # mass_spec = ContSpectra(mz=mass_x, intensities=mass_spec, metadata={'normalized': normalize})
    # add RT ranges to mass_spec metadata
    # mass_spec.set('rt_ranges', np.array([ranges_start, ranges_end]).T)
    return mass_spec_df  # Attention, return value changed [intensity] -> [intensity_norm]


def get_RT_spec_for_mass_range(data, ranges_start:np.array = None, ranges_end:np.array = None, rt_index_to_rt_mode:int = 0, rt_index_to_rt_params:list = [float(1/5),0], mass_index_to_mass_mode: int = 0, mass_index_to_mass_params:list = [2835.23579263574, -4049.42781522328], range_tol: float=0.1, range_mode:str = 'abs', sum_x:str = 'mass', normalize:bool = False):
    """
    This function returns the RT spec for the given mass range{s}.
    Parameters:
    data: GC-MS data
    ranges_start: list of start values of the ranges (in mass or index units, depending on the range_mode)
    ranges_end: list of end values of the ranges (mass units or index units, depending on the range_mode)
    rt_index_to_rt_mode: mode of the conversion from RT index to RT value. 0: linear conversion
    rt_index_to_rt_params: list of parameters for the conversion from RT index to RT value. The first parameter is the number of seconds per bin, the second parameter is the offset.
    mass_index_to_mass_mode: mode of the conversion from mass index to mass value. 0: square root conversion
    mass_index_to_mass_params: list of parameters for the conversion from mass index to mass value. The first parameter is the number of seconds per bin, the second parameter is the offset.
    range_tol: tolerance for the range values (in mass or index units, depending on the range_mode)
    range_mode: mode of the range values. 'abs': absolute values, 'ind': index values
    sum_x: x axis of the sum. 'RT': RT values, 'RT_index': RT index values
    normalize: if True, normalize the RT spec
    Returns:
    RT_spec_df: pandas dataframe with the RT spec data. The columns are 'RT' and 'intensity' (or 'intensity_norm' if normalize is True) 
    """
    if range_mode == 'abs' and ranges_start is not None and ranges_end is not None:
        # convert masses to mass indici
        ranges_start, _ = mass_to_mass_index(ranges_start, params=mass_index_to_mass_params, mode=mass_index_to_mass_mode)
        ranges_end, _ = mass_to_mass_index(ranges_end, params=mass_index_to_mass_params, mode=mass_index_to_mass_mode)
        range_tol, _ = mass_to_mass_index(range_tol, params=mass_index_to_mass_params, mode=mass_index_to_mass_mode)
        # convert to int
        ranges_start = ranges_start.astype(int)
        ranges_end = ranges_end.astype(int)
        # range_tol = range_tol.astype(int)
    # get the sum of the data in the given range
    RT_spec = get_sum(data, sum_type="RT", range_starts=ranges_start, range_ends=ranges_end, tol=range_tol)
    # convert the RT indici to RT values
    rt_index = np.arange(RT_spec.shape[0])
    if sum_x == 'RT':
        rt_x, _ = rt_index_to_rt(rt_index, params=rt_index_to_rt_params, mode=rt_index_to_rt_mode)
    elif sum_x == 'RT_index':
        rt_x = rt_index
    else:
        print("Error: Invalid sum_x. Choose from 'RT' or 'RT_index'")
        return None
    #if normalize:
    #    RT_spec = RT_spec/RT_spec.sum()
    # convert to pandas dataframe
    # generate data
    RT_spec_df = pd.DataFrame()
    RT_spec_df['RT'] = rt_x
    RT_spec_df['intensity'] = RT_spec
    #data_RT_spec = np.array([rt_x, RT_spec]).T
    #RT_spec_df = pd.DataFrame(data_RT_spec, columns=[sum_x, 'intensity'])
    if normalize:
        #RT_spec_df = RT_spec_df.rename(columns={'intensity': 'intensity_norm'})
        RT_spec_df['intensity_norm'] = RT_spec/RT_spec.sum()
    # convert to ContSpec object # TODO: write such a class!
    # RT_spec = ContSpectra(mz=rt_x, intensities=RT_spec, metadata={'normalized': normalize})
    # add mass ranges to RT_spec metadata
    #RT_spec.set('mass_ranges', np.array([ranges_start, ranges_end]).T)
    return RT_spec_df  # Attention, return value changed [intensity] -> [intensity_norm]

def get_ranges_for_summing(data, sum_type="RT", threshold_median_multiples=2, thresh_min_range_len = 3):
    """ This function returns the ranges for summing the data in the given sum_type by returning all ranges above a certain treshold.
    Parameters:
    data: GC-MS data to be summed
    sum_type: "RT" (get the mass index ranges, with have highest contribution to RT sum spec) or "mass" (get the RT index ranges which have the highest contribution to the summed mass spec)
    threshold_median_multiples: threshold for the sum of the data. The threshold is calculated as the median of the sum of the data multiplied by this value.
    thresh_min_range_len: minimum length of the ranges
    Returns:
    range_starts: list of start values of the ranges in the mass regime (if sum_type is "RT") or in the RT regime (if sum_type is "mass")
    range_ends: list of end values of the ranges in the mass regime (if sum_type is "RT") or in the RT regime (if sum_type is "mass")
    """
    if sum_type == "RT":
        summed_mass_val_median = np.median(data.sum(axis=1)) #sum over all RTs
        # arange values for RT
        int_values = np.arange(data.shape[1])
        # get mask for int values that are above the threshold, i.e. where the summed mass values are above the threshold
        range_mask = data.sum(axis=0) >= threshold_median_multiples * summed_mass_val_median

        range_starts = []
        range_ends = []
        range_start = None
        for i, mask_value in enumerate(range_mask):
            if mask_value and range_start is None:
                range_start = i
            elif not mask_value and range_start is not None:
                range_end = i-1
                range_starts.append(range_start)
                range_ends.append(range_end)
                range_start = None
        if range_start is not None:
            range_starts.append(range_start)
            range_ends.append(data.shape[1]-1)
        range_starts = np.array(range_starts)
        range_ends = np.array(range_ends)
        ind_above_thresh = np.where(range_ends-range_starts > thresh_min_range_len)
        range_starts = range_starts[ind_above_thresh]
        range_ends = range_ends[ind_above_thresh]
        return range_starts, range_ends
    elif sum_type == "mass":
        summed_rt_val_median = np.median(data.sum(axis=0)) #sum over all masses
        # arange values for mass
        int_values = np.arange(data.shape[0])
        # get mask for int values that are above the threshold, i.e. where the summed RT values are above the threshold
        range_mask = data.sum(axis=1) >= threshold_median_multiples * summed_rt_val_median
        # get the ranges
        range_starts = []
        range_ends = []
        range_start = None
        for i, mask_value in enumerate(range_mask):
            if mask_value and range_start is None:
                range_start = i
            elif not mask_value and range_start is not None:
                range_end = i-1
                range_starts.append(range_start)
                range_ends.append(range_end)
                range_start = None
        if range_start is not None:
            range_starts.append(range_start)
            range_ends.append(data.shape[0]-1)
        range_starts = np.array(range_starts)
        range_ends = np.array(range_ends)
        ind_above_thresh = np.where(range_ends-range_starts > thresh_min_range_len)
        range_starts = range_starts[ind_above_thresh]
        range_ends = range_ends[ind_above_thresh]
        return range_starts, range_ends
    else:
        print("Error: Invalid sum type. Choose from 'RT' or 'mass'")
        return None

def combinations_sum(target, candidates):
    """
    Given a set of candidate numbers (candidates) (without duplicates) and a target number (target),
    find all unique combinations in candidates where the candidate numbers sum to target.
    The same repeated number may be chosen from candidates unlimited number of times.
    Note:
    All numbers (including target) will be positive integers.
    The solution set must not contain duplicate combinations.
    Parameters:
    target: target number
    candidates: list of candidate numbers
    Returns:
    result: list of lists of the unique combinations, where each inner list is a combination that sums to target
    """
    def backtrack(start, target, path):
        if target == 0:
            result.append(path)
            return
        for i in range(start, len(candidates)):
            if candidates[i] > target:
                break
            backtrack(i, target - candidates[i], path + [candidates[i]])

    result = []
    candidates.sort()
    backtrack(0, target, [])
    return result

def get_unit_mass_of_largest_isotope(elem):
    """
    Get the mass of the highest isotope of an element.
    Parameters:
    elem: element as string
    Returns:
    mass: unit mass of the highest isotope of the element
    """
    isospec = chem_formula(elem).spectrum()
    return [int(mz) for mz in isospec if isclose(isospec[mz].intensity, 100, rel_tol = 1E-9)][0]
# get mass of highest isotope (rel equivalent to 100%)(just to be sure, a bit complicated)

def check_dbe_senior_rule_and_if_existence_is_likely(chem_formula_i:chem_formula) -> list[bool, float]:
    """ See: alpinac.periodic_table: test the SENIOR theorem, item (ii),
        as described in Kind and Fiehn 2007:
        The sum of valences is greater than or equal to twice the
        maximum valence.
        This prevents fragments such as CFCl to be considered as
        valid molecular ion.
        THe DBE value is connected to the ability of a molecule to bind to more H atoms (open bound sites).
        The DBE value is calculated as follows:
        DBE = 2 + (sum of valences - 2*max valence)/2
        Parameters:
        chem_formula_i: chem_formula object
        Returns:
        dbe: double bond equivalent
        Senior_rule_test: True if the Senior rule is fulfilled
        likely: True if the Senior rule is fulfilled and the existence of the molecule is likely by dbe.
        """
    composition =  chem_formula_i.composition()
    sum_valence = 0
    max_valence = 0
    dbe_part = 2

    for line in composition:
        atom_count = composition[line].count
        atom_only_largest_isotope = ''.join([i for i in line if i.isalpha()])  # valence can not be calculated for [13]C, ...
        atom_valence = abs(min(element(atom_only_largest_isotope).oxistates)) # ground-state valence is absolute of lowest state
        dbe_part += atom_count * (atom_valence - 2)
        sum_valence += atom_count * (atom_valence)
        max_valence = max(max_valence, atom_valence)
    dbe = dbe_part/2
    Senior_rule_test = sum_valence >= 2*max_valence
    likely = Senior_rule_test and dbe >= 0
    return dbe, Senior_rule_test, likely

def get_chemical_possible_combinations(elements: List[str], target_masses: int or list, radicals_allowed: bool = False) -> List[chem_formula]:
    """ Get all possible combinations of unit masses that sum up to target_masses (one specific or a range).
    Parameters:
        elements: list of elements as strings
        target_masses: list of target masses as integers or a single target mass as integer.
        radicals_allowed: if True, radicals are allowed (for example as fragment, use for molecular ion in EI: False )
        Returns:
        formulas_res: list of formulas in ground-state
        formulas_rad: list of potential radical formulas
    """
    #deep copy of elements
    #elements_tmp = elements.copy()
    if isinstance(target_masses, int):
        target_masses = [target_masses]

    formulas_res = []
    formulas_rad = []
    unit_masses = [get_unit_mass_of_largest_isotope(elem) for elem in elements]
    dict_unit_mass_and_element = {key: value for key, value in zip(unit_masses, elements)}
    for target_mass in target_masses:
        elements_tmp = np.array(elements).tolist() # deep copy of elements
        #flag = 1
        res = combinations_sum(target = target_mass, candidates = unit_masses)
        if len(res) > 0:
            for row in res:
                # test if the combination is possible
                # count the number of atoms of each element
                occurences_numbers = Counter(row)
                occurences_atoms = {dict_unit_mass_and_element[key]: value for key, value in occurences_numbers.items()}
                # get the element of each unit mass
                formula_parts = [str(atom) + str(occurences_atoms[atom])  for atom in occurences_atoms.keys()]
                formula_row = chem_formula("".join(formula_parts))
                # if can exist chemically:
                if check_dbe_senior_rule_and_if_existence_is_likely(formula_row)[2]:
                    #if flag == 1:
                    #    print("target_mass: {}".format(target_mass))
                    #    flag = 0                   
                    print(formula_row.formula)
                    formulas_res.append(formula_row.formula)
                # TODO check octed rules for IONS! Molecules where the ion gets more stable by adding a H atom, or Cl,F,... because of the octed rule, are more likely to exist.
                # so if molecule is formed which would have a closed shell due to one elctron missing, it is more likely to exist.
                if radicals_allowed: # if radicals are allowed this correspond to a stable formula + H atom (just that the H atom is not there, but a free bound
                    # -> we test if the formula + H atom is stable: this would be a radical then
                    # as long as DBE is positive/zero H atoms can be added # check this!
                    # if 'H' in elements add H atom
                    dbe = 0
                    while dbe >= 0:
                        if 'H' in occurences_atoms:
                            occurences_atoms['H'] += 1
                            radical_parts = [str(atom) + str(occurences_atoms[atom])  for atom in occurences_atoms.keys()]
                            radical_row = chem_formula("".join(radical_parts))
                        else:
                            # add 'H' to elements
                            elements_tmp.append('H')
                            occurences_atoms['H'] = 1
                            radical_parts = [str(atom) + str(occurences_atoms[atom])  for atom in occurences_atoms.keys()]
                            radical_row = chem_formula("".join(radical_parts))
                        # generate fake formula
                        if check_dbe_senior_rule_and_if_existence_is_likely(radical_row)[2]:
                        #if flag == 1:
                        #    print("target_mass: {}".format(target_mass))
                        #    flag = 0
                            print('radical:')
                            print(formula_row.formula)
                            formulas_rad.append(formula_row.formula)
                            break
                        dbe = check_dbe_senior_rule_and_if_existence_is_likely(radical_row)[0]     
        # delete elements_tmp, as not needed anymore
        del elements_tmp                         
    return np.array(formulas_res) , np.array(formulas_rad)
    
def get_calibrant_peaks(data, start_or_end="start", factor_median_for_threshold=10, plot=False):
    """ This function returns the index of the calibrant peak in the data.
    Parameters:
    data: data as numpy array
    start_or_end: "start" or "end". If "start", the first calibrant peak set is returned. If "end", the last calibrant peak set is returned.
    tol: tolerance for the range values
    nbr_cal_peaks: number of calibrant peaks to be returned
    returns:
    left_rising_edge_index: index of the left rising edge of the calibrant peak
    right_falling_edge_index: index of the right falling edge of the calibrant peak
    peaks: index of the calibrant peak
    rt_sum_threshold_index_interval_index[flank_indici]: index of the flank of the calibrant peak
    rt_sum_threshold_index_interval_index[half_max_index_rising]: index of the left half maximum of the calibrant peak
    rt_sum_threshold_index_interval_index[half_max_index_falling]: index of the right half maximum of the calibrant peak
    rt_index_to_rt_flank: RT of the flank of the calibrant peak
    rt_index_to_rt_leftFWHM: RT of the left half maximum of the calibrant peak
    rt_index_to_rt_rightFWHM:RT of the right half maximum of the calibrant peak
    """
    if start_or_end == "start":
        ind = 0
        low_range_for_find_peaks = -20
        high_range_for_find_peaks = 150
    elif start_or_end == "end":
        ind = -1
        low_range_for_find_peaks = -150
        high_range_for_find_peaks = 20


    # get calibrant peaks
    # search for first maxima in retention time sum plot
    rt_sum = get_sum(data, sum_type="RT")
    #rt_sum = rt_sum.sum(axis=1)
    # find first time where the sum is larger than 10 times the median
    rt_sum_median = np.median(rt_sum)
    rt_sum_threshold = factor_median_for_threshold* rt_sum_median
    rt_sum_threshold_indices = np.where(rt_sum > rt_sum_threshold)[ind]
    rt_sum_threshold_index = rt_sum_threshold_indices[ind]
    # get data in an interval of -1 +10s around that index
    rt_sum_threshold_index_start = rt_sum_threshold_index + low_range_for_find_peaks
    rt_sum_threshold_index_end = rt_sum_threshold_index + high_range_for_find_peaks
    rt_sum_threshold_index_start = max(0, rt_sum_threshold_index_start)
    rt_sum_threshold_index_end = min(rt_sum_threshold_index_end, rt_sum.shape[0])
    # plot rt_sum in that interval
    rt_sum_threshold_index_interval = rt_sum[rt_sum_threshold_index_start:rt_sum_threshold_index_end]
    rt_sum_threshold_index_interval_index = np.arange(rt_sum_threshold_index_start, rt_sum_threshold_index_end)
    
    # find peaks in that interval
    peaks, _ = find_peaks(rt_sum_threshold_index_interval, height=0.15*np.max(rt_sum_threshold_index_interval), prominence=0.15*np.max(rt_sum_threshold_index_interval))
    if plot:
        fig, ax = plt.subplots(1,1)
        plt.plot(rt_index_to_rt(rt_sum_threshold_index_interval_index)[0], rt_sum_threshold_index_interval)
        # add to plot
        plt.plot(rt_index_to_rt(rt_sum_threshold_index_interval_index[peaks])[0], rt_sum_threshold_index_interval[peaks], "x")
    # seperate rt_sum_threshold_index_interval into 4 windows by detecting where a flank is rising
    flank_indici = []
    left_rising_edge_index = []
    right_falling_edge_index = []
    # go from each maximum down until signal is rising again
    for i in range(len(peaks)):
        # get the index of the first peak
        peak_index = peaks[i]
        # get the index of the first flank
        flank_index = np.where(rt_sum_threshold_index_interval[:peak_index] < rt_sum_threshold_index_interval[peak_index])[0][-1]
        flank_indici.append(flank_index)
        # plot the point where the flank is rising as a vertical line
        # get the index where the signal is fallen to 50% of the maximum
        half_max = 0.5*rt_sum_threshold_index_interval[peak_index]
        # get the first value when signal fall to half_max after the peak
        half_max_index_falling = np.where(rt_sum_threshold_index_interval[peak_index:] < half_max)[0][0] + peak_index
        # get the first value when signal fall to half_max before the peak
        half_max_index_rising = np.where(rt_sum_threshold_index_interval[:peak_index] < half_max)[0][-1]
        #append to left_rising_edge_index
        left_rising_edge_index.append(half_max_index_rising)    
        #append to right_falling_edge_index
        right_falling_edge_index.append(half_max_index_falling)

        if plot:
            # plot the point where the signal is at 50% of the maximum as a vertical line
            plt.axvline(rt_index_to_rt(peaks[i] + rt_sum_threshold_index_start)[0], color="red", linestyle="--", linewidth=0.5)
            plt.axvline(rt_index_to_rt(half_max_index_rising + rt_sum_threshold_index_start)[0], color="green", linestyle="--", linewidth=0.5)
            plt.axvline(rt_index_to_rt(half_max_index_falling + rt_sum_threshold_index_start)[0], color="green", linestyle="--", linewidth=0.5)
        # append to rt_sum_threshold_index_interval_indici
    # correct for the offset of the interval
    left_rising_edge_index = np.array(left_rising_edge_index) + rt_sum_threshold_index_start
    right_falling_edge_index = np.array(right_falling_edge_index) + rt_sum_threshold_index_start
    peaks = np.array(peaks) + rt_sum_threshold_index_start

    
    return_list = [left_rising_edge_index,
                    right_falling_edge_index,
                    peaks, 
                    rt_index_to_rt(left_rising_edge_index)[0],
                    rt_index_to_rt(right_falling_edge_index)[0],
                    rt_index_to_rt(peaks)[0]]
    return return_list

    # TODO: check if target mass already in the range of masses and extend if not

def get_chemical_possible_combinations_for_elem_list(elements, target_masses: np.array, radicals_allowed:bool = False, path_to_files:str = None, recalculate:bool = False):
    """ Get all possible combinations of unit masses that sum up to target_masses (one specific or a range).
        :par: elements: list of elements as strings
        :par: target_masses: list of target masses as integers or a single target mass as integer.
        :return: formulas_gs: list of formulas in ground-state
        :return: formulas_radicals: list of potential radical formulas
        :return: formulas: list of all formulas
        :return: masses: list of masses
    """
    # deep copy of elements
    elements_tmp = elements.copy()
    #target_mass_str = [str(target_mass) for target_mass in target_masses]
    masses = np.array([])
    formulas = np.array([])
    masses_no_combi = np.array([])

    # check if file formulas_elements.txt exists
    if path_to_files is None:
        path_to_file = str("formulas_"+"".join(elements)+"_radicals_" + str(radicals_allowed) + ".txt")
    else:
        path_to_file = path_to_files + "/" + str("formulas_"+"".join(elements)+"_radicals_" + str(radicals_allowed) + ".txt")

    if recalculate:
        #delete path_to_file if it exists
        if Path(path_to_file).exists():
            os.remove(path_to_file)

    if Path(path_to_file).exists() and not recalculate:
        data = pd.read_csv(path_to_file, sep="\t", names = ["formula", "mass"])
        formulas = data["formula"].values
        masses = data["mass"].values
        # get masses for which previous calculations did not yield a result:
        path_to_file_masses_with_no_results = path_to_file[:-4] + "_masses_with_no_combi.txt"
        if Path(path_to_file_masses_with_no_results).exists():
            masses_no_combi = pd.read_csv(path_to_file_masses_with_no_results, sep="\t", names = ["mass"])["mass"].values
        else:
            masses_no_combi = np.array([])
        # opposite of intersect between masses and target_masses
        masses_not_in_masses = np.setdiff1d(target_masses, masses.astype(int))
        masses_to_calculate = np.setdiff1d(masses_not_in_masses, masses_no_combi)       


        formulas_gs = []
        formulas_radicals = []
        # determine 
    else:
        masses_to_calculate = target_masses
    if len(masses_to_calculate) > 0:
        print("calculate combinations for masses: {}".format(masses_to_calculate))
        # calculate with backtrack all possible combinations of unit masses that sum up to each element of the target masses
        formulas_gs, formulas_radicals = get_chemical_possible_combinations(elements = elements_tmp, target_masses = masses_to_calculate, radicals_allowed = radicals_allowed)
        # join ndarrays
        formulas_calc = np.concatenate((formulas_gs, formulas_radicals))
        # get unique formulas
        formulas_calc = np.unique(formulas_calc)    
        # get masses as np array
        masses_calc = np.array([chem_formula(formula).isotope.mass for formula in formulas_calc])
        # sort masses and formulas by mass
        #ind_sort = np.argsort(masses_calc)
        #masses_calc = masses_calc[ind_sort]
        #formulas_calc = formulas_calc[ind_sort]
        # add to ndaraay masses and ndarray formulas
        masses = np.concatenate((masses, masses_calc))
        formulas = np.concatenate((formulas, formulas_calc))
        # sort masses and formulas by masses:
        ind_sort = np.argsort(masses)
        masses = masses[ind_sort]
        formulas = formulas[ind_sort]   


        #formulas_calc.sort(key=lambda x: chem_formula(x).isotope.mass)
        # get sorted masses
        #masses_calc = [round(chem_formula(formula).isotope.mass, ndigits=8) for formula in formulas_calc]
        # add to masses and formulas
        #masses.extend(masses_calc)
        #formulas.extend(formulas_calc)
        if len(masses) > 0:
            # to pandas dataframe with precision 8
            data_updated = pd.DataFrame({"formula": formulas, "mass": masses})
            # get masses with no results
            masses_not_calculable = np.setdiff1d(target_masses, masses.astype(int))
            # add to masses_no_combi
            masses_no_combi = np.concatenate((masses_no_combi, masses_not_calculable))
            # write_it_to_file
            if len(masses_no_combi) > 0:
                path_to_file_masses_with_no_results = path_to_file[:-4] + "_masses_with_no_combi.txt"
                pd.DataFrame(masses_no_combi).to_csv(path_to_file_masses_with_no_results, sep="\t", index=False, header=False)

            # save as csv
            data_updated.to_csv(path_to_file, sep="\t", index=False, header=False)
    return formulas_gs, formulas_radicals, formulas, masses

def calculate_interpolated_intensities(masses, formulas, mass_indici, intensities, threshold = 0):
    """
    Calculate the intensities of the mass spectra at the masses of the formulas.
    Parameters:
    masses: list of masses
    formulas: list of formulas
    mass_sum_all_indici: array of indices of the mass spectra
    mass_sum_all: array of mass spectra 
    """

    # for each intensity_list_elem calculate the intensities
        # interpolate the intensity at the masses of the formulas
    intensities_interp = np.interp(x = masses, xp = index_to_mass(mass_indici), fp = intensities)
    mask = intensities_interp >= threshold
    fig, ax = plt.subplots(1,1)
    plt.plot( index_to_mass(mass_indici), intensities)
    plt.plot(masses, intensities_interp, "x")
    plt.show()

    return intensities_interp[mask], np.array(formulas)[mask], np.array(masses)[mask]

def get_mass(spikes: list) -> (list, list):
    """ Get the mass of a spike
    Parameters:
    spikes: list of spikes
    """
    masses = []
    formulas = []
    for i, spike in enumerate(spikes):
        chem_form = chem_formula(spike)
        formulas.append(chem_form)
        masses.append(chem_formula(spike).isotope.mass)
    return formulas, masses

def index_to_mass(i, p1, p2):
    """
    Convert index to mass
    Parameters:
    i: index
    p1: parameter 1
    p2: parameter 2
    """
    return ((np.array(i) - p2) / p1) ** 2

def calculate_interpolated_intensities(masses, formulas, mass_indici, intensities, params_initial, threshold=0, index_to_mass = index_to_mass, plot=False):
    """
    Calculate the intensities of the mass spectra at the masses of the formulas.
    Parameters:
    masses: list of masses
    formulas: list of formulas
    mass_sum_all_indici: array of indices of the mass spectra
    mass_sum_all: array of mass spectra
    params_initial: initial parameters for the optimization
    threshold: threshold for the intensities
    index_to_mass: function to convert index to mass
    plot: if True, plot the results
    """
    x0 = params_initial

    def objective(p):
        # Calculate interpolated intensities
        intensities_interp = np.interp(x=masses, xp=index_to_mass(mass_indici, *p), fp=intensities)
        # Find indices of intensities above threshold
        mask = intensities_interp > threshold
        # Return negative number of intensities and negative total intensity
        # return -(mask.sum() + intensities_interp[mask].sum()/intensities_interp.sum())
        # return ( intensities_interp[mask].sum())
        return (mask.sum()/len(masses) + intensities_interp[mask].sum()/intensities.sum())


    def callback(xk):
        # Check if new parameters fall within allowed variation from initial parameters
        if not all((np.abs(xk - x0) / np.abs(x0)) < 0.1):
            print("Warning: parameters outside of variation limit")


    # Minimize negative objective to maximize number of intensities and total intensity
    res = minimize(lambda x: -objective(x), x0, method='Nelder-Mead', callback = callback, bounds = bounds)
    print(res)

    # Calculate interpolated intensities with optimized parameters
    p_opt = res.x
    intensities_interp = np.interp(x=index_to_mass(mass_indici, *p_opt), xp=index_to_mass(mass_indici, *p_opt), fp=intensities)
    # masses optimized
    masses_opt = index_to_mass(mass_indici, *p_opt)
    # plot masses optimized vs intensities
    fig, ax = plt.subplots(1, 1)
    plt.plot(index_to_mass(mass_indici, *x0), intensities)
    plt.plot(masses_opt, intensities_interp, "red")
    # Find indices of intensities above threshold
    for mass in masses:
        plt.axvline(x = mass,  color="red", linestyle="--", linewidth=0.5)


    mask = intensities_interp >= threshold
    return intensities_interp[mask], np.array(masses_opt)[mask], res


def clean_empa_name_string(name):
    """Clean the name string of the empa documents to get the formula string
    Parameters:
    name: name string
    Returns:
    formula_str: formula string
    """
    formula_str = name.replace("_myriam", "").replace("_nist", "").replace("_xtra_CH", "").strip()
    formula_str = formula_str.replace(" (e)", "").replace(" (n)", "").replace("_s", "").replace("xtra_CH", "").replace("-only", "").strip()
    # if ends with "E" or "Z" replace it with "(E)-" or "(Z)-" at the beginning of the string
    if formula_str.endswith("E") or formula_str.endswith("Z"):
        formula_str = "(" + formula_str[-1] + ")-" + formula_str[:-1]
    # if formula starts with "c-" replace it with "cyclo"
    if formula_str.startswith("c-"):
        formula_str = "cyclo" + formula_str[2:]
    # if formula starts with "H-" followed by four numbers replace it with "Halon-" followed by the four numbers
    if formula_str.startswith("H-") and len(formula_str) == 6:
        formula_str = "Halon-" + formula_str[2:]
    return formula_str


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


def formula_by_name(name:str, formula_dict:dict, formulas_dict_manually_added:dict = {}, target_list = []):
    """
    This function returns the formula of a compound by its name
    :param name: name of the compound (name as saved in match-library)
    :return: formula of the compound
    :return: iupac_name of the compound
    :return: cid of the compound
    :return: inchi of the compound
    :return: RT_expected of the compound
    """
    #formula_dict
    iupac_name = name # if we have nothing better, we use the name as iupac_name
    cid = np.nan
    inchi = np.nan
    RT_expected = np.nan
    formula_str = clean_empa_name_string(name)
    numerical_formula = False
    #check if string is already a key of the dict formula_dict
    if formula_str in formula_dict:
        formula, iupac_name , cid, inchi, RT_expected = formula_dict[formula_str]
    else:
        try: # check if formula_str is a numerical formula
            formula = chem_formula(formula_str).formula
            numerical_formula = True
        except:
            numerical_formula = False
        # check if in manual added formulas
        if formula_str in formulas_dict_manually_added:
            formula, iupac_name ,cid, inchi, RT_expected = formulas_dict_manually_added[formula_str]
        # check if isinstance chem_formula, if yes formula_dict = str(chem_formula(formula_str)).formula
        elif numerical_formula:
            formula = str(chem_formula(formula_str).formula)
            try:
                results = pcp.get_compounds(formula, 'formula')
                if len(results) == 1: # one formula has normally many results unless very simple (Xe, ...)
                    compound = results[0]
                    iupac_name = compound.iupac_name
                    cid = str(compound.cid)
                    inchi = str(compound.inchi)
                    if inchi in target_list['inchi_key'].keys():
                        RT_expected = target_list.loc[target_list['inchi_key'] == inchi, 'RT'].values[0]
                    # add to dict if not already in it
            except:
                print('no unique compound found for formula: {}'.format(formula))
            formula_dict[formula_str] = [formula, iupac_name, str(cid), str(inchi), str(RT_expected)]
        # check if cas number
        else:
            try:
                if formula_str.replace("-", "").replace(" ", "").isnumeric():
            # get inchi by cas number using cirpy
                    inchi = cirpy.resolve(formula_str, 'inchi')
                    results = pcp.get_compounds(inchi, 'inchi')
                else:
                    # Search for compounds with the name, e.g. "ethane"
                    results = pcp.get_compounds(formula_str, 'name')
                # Get the first result from the search
                compound = results[0]
                # Get the molecular formula from the compound
                formula = compound.molecular_formula
                iupac_name = str(compound.iupac_name)
                cid = str(compound.cid)
                inchi = str(compound.inchi)

            except:
                try:
                    formula = halocarbname2formula(formula_str)[0]
                except:
                    formula = ""
        
            if len(target_list)>0 and not isinstance(inchi, str):
                if inchi in target_list['inchi_key'].values:
                    print("found in target list")
                    RT_expected = target_list.loc[target_list['inchi_key'] == inchi, 'RT'].values[0]
            # add to dict if not already in it
            formula_dict[formula_str] = [formula, iupac_name, str(cid), str(inchi), str(RT_expected)]
            # check if RT is available in Target List by comparing inchi_keys
    
    return formula, iupac_name, cid, inchi, RT_expected, formula_dict




# excecute if called as main
if __name__ == "__main__":
    # all integer masses from 1 to 300
    masses = np.arange(1, 300)
    #get_chemical_possible_combinations_for_elem_list(['C', 'Cl', 'F'], masses, radicals_allowed=False, path_to_files=r"C:\Users\kaho\polybox\Measurement_Campaign\formulas_non_target", recalculate=True)
    get_chemical_possible_combinations_for_elem_list(['C', 'Cl', 'F', 'Br', 'H', 'N', 'O', "I", 'S', 'P'], masses, radicals_allowed=False, path_to_files=r"C:\Users\kaho\polybox\Measurement_Campaign\formulas_non_target", recalculate=True)

    stop()


    # all integer masses from 1 to 300
    masses = np.arange(1, 300)
    #get_chemical_possible_combinations_for_elem_list(['C', 'Cl', 'F'], masses, radicals_allowed=False, path_to_files=r"C:\Users\kaho\polybox\Measurement_Campaign\formulas_non_target", recalculate=True)
    #get_chemical_possible_combinations_for_elem_list(['C', 'Cl', 'F', 'Br', 'H'], masses, radicals_allowed=False, path_to_files=r"C:\Users\kaho\polybox\Measurement_Campaign\formulas_non_target", recalculate=True)
    get_chemical_possible_combinations_for_elem_list(['C', 'Cl', 'F', 'Br', 'H', 'N', 'O', 'I', 'P', 'S'], masses, radicals_allowed=False, path_to_files=r"C:\Users\lab134\Desktop\ecToFKampagne\Gemessene_Proben", recalculate=True)

    stop()

    formula_by_name('463-49-0',{})
    # LOAD DATA
    file = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230401.1747.air.3.h5")
    file = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303_interesting_file\230311.0125.tank.1.h5")
    ei_seg = (3,4)
    ci_seg = (1,2)
    ionisation_dict = {"EI": ei_seg,"CI": ci_seg}
    # These are some masses that we don't want as they have maybe a column bleed
    masses_to_ignore = [185, 187]

    data_EI, reader_EI = load_data(file, ionisation_dict, 'EI')
    sumd = get_sum(data_EI, sum_type = 'mass', tol=0.1)

    stop()

    data_CI, reader_CI = load_data(file, ionisation_dict, 'CI')

    
    # get mass cailbration from hdf5 file
    # manual read out of params
    # Myriams calibration
    #reader.mass_calibration_data['EI'].mass_cal_parameters
    # TOFWerk calibration
    #reader.mass_calibration_data['EI'].get_parameters_for_rt(125)
    #reader.mass_cal_tofwerk
    #reader.rt_calibration_data
    


    # get mass calibration from hdf5 file
    #read_mass_cal_par(reader)
    # i = p1*np.sqrt(mz) + p2
    p1 = 2835.23579263574
    p2 = -4049.42781522328

    #RT calibration
    s_per_bin = float(1/5) # 5Hz
    offset = 0 

    

    data = data_EI
    ranges_start_example = [600]
    ranges_end_example = [700]
    mass_spec_df= get_mass_spec_for_RT_range(data, ranges_start = ranges_start_example, ranges_end = ranges_end_example, rt_index_to_rt_mode = 0, rt_index_to_rt_params=[float(1/5),0], mass_index_to_mass_mode = 0, mass_index_to_mass_params=[2835.23579263574, -4049.42781522328], range_tol=0.1, range_mode='abs', sum_x='mass', normalize=False) 

    fig, ax = plt.subplots(1,1)
    ax.plot(mass_spec_df['mass'], mass_spec_df['intensity'])
    #add metadata as title
    plt.title("RT range: {} - {}".format(ranges_start_example, ranges_end_example   ))


    ranges_start_example = [20]
    ranges_end_example = [150]
    RT_spec_df = get_RT_spec_for_mass_range(data, ranges_start = None, ranges_end = None, rt_index_to_rt_mode = 0, rt_index_to_rt_params=[float(1/5),0], mass_index_to_mass_mode = 0, mass_index_to_mass_params=[2835.23579263574, -4049.42781522328], range_tol=0.1, range_mode='abs', sum_x='RT', normalize=False)
    fig, ax = plt.subplots(1,1)
    ax.plot(RT_spec_df['RT'], RT_spec_df['intensity'])
    #add metadata as title
    plt.title("mass range: {} - {}".format(ranges_start_example, ranges_end_example))


    # add compounds vertically written to the plot
    #read_determined_peaks_data
    path_det_peak_data = Path(r"c:\Users\kaho\polybox\Presentations\PythonPlots\peaks_det_230211.0125.tank.1.xlsx")
    det_peak_data = pd.read_excel(path_det_peak_data)
    # delete all rows which are Nan in first column
    det_peak_data = det_peak_data[det_peak_data.iloc[:,0].notna()]
    # delete rows where first column is "no usable data"
    det_peak_data = det_peak_data[det_peak_data.iloc[:,0] != "no usable data"]
    # split second column by ;\n in two columns
    det_peak_data[['RT_s', 'RT_min']] = det_peak_data.iloc[:,1].str.split(';\n',  expand=True)
    pos_label = det_peak_data['RT_s']
    # convert to float
    pos_label = pos_label.astype(float)
    label_name = det_peak_data['eval_status'].astype(str)
    # convert to np.array
    pos_label = np.array(pos_label)
    label_name = np.array(label_name)
    # add to plot as vertical labels at 0
    for i, label in enumerate(label_name):
        ax.text(pos_label[i], 0, label, rotation=90, color='red', fontsize=8, verticalalignment='top', horizontalalignment='center')

    # get maximum intensity




    # get first 4 and last 4 maxima (spikes)
    # this can be extracted from TPSscripting data and comparison to the peaks found in peakfinder
    # the calibrant trigger is at:
    calibrant_theo = [120,125,130,135,3065,3070,3075,3080]
    # number of calibrant peaks can be extracted as well

    p1 = 2835.23579263574
    p2 = -4049.42781522328


    # for fit:
    # Initial guess for p1 and p2
    p1_init = 2835.23579263574
    p2_init = -4049.42781522328
    x0 = np.array([p1_init, p2_init])
    # Define parameter bounds
    p1_bounds = (p1_init - 0.002*p1_init, p1_init + 0.002*p1_init)
    p2_bounds = (p2_init - 0.002*p2_init, p2_init + 0.002*p2_init)
    # sort bounds
    if p1_bounds[0] > p1_bounds[1]:
        p1_bounds = (p1_bounds[1], p1_bounds[0])
    if p2_bounds[0] > p2_bounds[1]:
        p2_bounds = (p2_bounds[1], p2_bounds[0])
    bounds = [p1_bounds, p2_bounds]

    # Define initial guess
    params_initial = np.array([(p1_bounds[0] + p1_bounds[1])/2, (p2_bounds[0] + p2_bounds[1])/2])



    data_CI_sum = get_sum(data_CI, "RT")
    # plot vs mass calibrated data

    
    # get sum:
    data_sum = get_sum(data, "mass")
    # plot vs mass calibrated data:
    fig, ax = plt.subplots(1,1)
    plt.plot(index_to_mass(range(len(data_sum)), p1=p1, p2 = p2), data_sum)
    rt_sum = get_sum(data, sum_type="RT", range_starts=[0, 3000], range_ends=[200, 400], tol=0.1)
    mass_sum = get_sum(data, sum_type="mass", range_starts=[10, 30], range_ends=[20, 40], tol=0.1)

    # sum all rt_sum data
    rt_sum_all = rt_sum.sum(axis=1)
    plt.plot(rt_sum_all)

    #fig, ax = plt.subplots(1,1)
    #plt.plot(rt_sum.sum)

    # get ranges for RT by choosing only the highest mass peaks


    rt_sum = get_sum(data, sum_type="RT", range_starts=[0], range_ends=[22248], tol=5)
    mass_sum = get_sum(data, sum_type="mass", tol=0.1)

    fig, ax = plt.subplots(1,1)
    plt.plot(mass_sum)

    # PLOT RT SUM


    # sum all rt_sum data
    fig, ax = plt.subplots(1,1)
    #rt_sum_all = rt_sum.sum(axis=1)
    # RT = index*s_per_bin + offset
    index =  np.arange(rt_sum.shape[0])  
    plt.plot(index*s_per_bin, np.array(rt_sum))


    # PLOT HIGHEST MASS PEAKS

    range_starts, range_ends = get_ranges_for_summing(data, sum_type="RT", threshold_median_multiples=2, thresh_min_range_len= 3)
    print(np.array(range_ends)-np.array(range_starts))
    # eliminiate ranges that are too small
    # print mass spectra for each range
    fig, ax = plt.subplots(1,1)
    range_center = 0.5*(range_starts + range_ends)
    line_type = ["-", "--", "-.", ":"]*100


    for i, range_start in enumerate(range_starts):
        range_end = range_ends[i]
        rt_sum = get_sum(data, sum_type="RT", range_starts=[range_start], range_ends=[range_end], tol=0.1)
        #plt.plot(rt_sum, label="RT range: "+str(range_center[i]))
        plt.plot(rt_sum, label=round(index_to_mass(range_center[i], p1=p1, p2=p2)), linestyle = line_type[i])
    # plot horizontal labels
    # if color palette is repeating change the line type
    plt.legend(bbox_to_anchor=(0.05, 1), loc='upper left', borderaxespad=0.0, fontsize=6, ncol=3)
    plt.show()

    # ONE COULD IMPLEMENT A RT TIME SEPARATION HERE


    # draw mass plot with ranges
    fig, ax = plt.subplots(1,1)
    plt.plot(data.sum(axis=0))
    for i, range_start in enumerate(range_starts):
        range_end = range_ends[i]
        plt.axvline(range_start, color="red", linestyle="--", linewidth=0.5)
        plt.axvline(range_end, color="lightgreen", linestyle="--", linewidth=0.5)
        # plot horizontal line at zero
        plt.axhline(0, color="black", linewidth=0.5)

  

    print(np.array(range_ends)-np.array(range_starts))

    
    range_start, range_end = get_ranges_for_summing(data, sum_type="RT", threshold_median_multiples=5000, thresh_min_range_len= 3)
    rt_sum = get_sum(data, sum_type="RT", range_starts=[range_start], range_ends=[range_end], tol=0.1)
    get_calibrant_peaks(data, start_or_end="end", factor_median_for_threshold=10, plot=True)
    get_calibrant_peaks(data, start_or_end="start", factor_median_for_threshold=10, plot=True)


    fig, ax = plt.subplots(1,1)
    plt.plot(data.sum(axis=0)) # axis 0: mass, 1: RT
    fig, ax = plt.subplots(1,1)
    plt.plot(rt_sum)


    # CALIBRATION

    # use calibrant peak to do a rough calibration



    left_rising_edge_index, right_falling_edge_index, peaks, _, _, _ = get_calibrant_peaks(data, start_or_end="start", plot=True)

    # sum mass spectrum in each window defined by left_rising_edge_index and right_falling_edge_index
    mass_sum = get_sum(data, sum_type="mass", range_starts=left_rising_edge_index, range_ends=right_falling_edge_index, tol=0.000001)
    # plot the sum of the mass spectrum
    fig, ax = plt.subplots(1,1)
    # plot sum of all mass spectra
    mass_sum_all = mass_sum.sum(axis=0)
    plt.plot(mass_sum_all)

    # detect all mass peaks
    peaks = find_peaks(mass_sum_all, height=0.005*np.max(mass_sum_all), prominence=0.005*np.max(mass_sum_all))
    #add to plot
    plt.plot(peaks[0], peaks[1]['prominences'], "x")



    spikes=[
            "CF3", 
            "C2F3",  "C2F4", 
            "C3F4", 
            "C4F5", 
            "C5F5", 
            "CF2","C3F3","CF", "C2F5", 
            "C3F5",
            # "C5NF10",
            "C3F6", "C3F7", "C4F7",
            "C4F9", 
            "C5F9", 
            "C6F9", 
            "C7F11", 
            "C8F9", 
            #"C6NF12", "C7NF14", "C8FN14", "C8NF16"
        ]
    #convert formula to counts of atoms
    spikes_counts = []
    for spike in spikes:
        spikes_counts.append(chem_formula(spike).atoms) 
    # get the sum of all spikes
    spikes_sum = np.sum(spikes_counts, axis=0)
    formulas, masses = get_mass(spikes)


    #1) test dbe
    #2) test SENIOR theorem
    # test_SENIOR_rule = sum([valence[i]*fragment[i] for i in range(len(fragment))]) >= 2 * max([valence[i] for i in range(len(fragment)) if fragment[i]>0])


    print(chem_formula("C2F3").composition())
    # for each atom in the composition calculate count times elem.oxistates - 2 and sum

    chem_formula_i = chem_formula("C2F6")



    print(check_dbe_senior_rule_and_if_existence_is_likely(chem_formula('H7N2O2')))

    print(chem_formula("C2F3").spectrum())
   

    C = element('C')
    # maximum of the oxistates of the element are the valence 
    print(max(C.oxistates))



    # check if file formulas_elements.txt exists
    #elements = ['C','H']
    #target_masses = np.arange(30, 300)

    

    # calculate the possible isotopic spectra of air constituents
    air_constitutents = ["O2", "N2", "Ar", "CO2", "H2O", "Xe", "Kr", "Ne", "CH4", "SO2", "NO", "NO2", "N2O", "NH3", "HNO3", "H2SO4", "HCl", "HBr", "HF", "H2S", "H2Se", "H2Te", "CO", "O3", "H2"]
    masses_air_constituents = []
    formulas_air_constituents = []
    for elem in air_constitutents:
        # get isotope spectrum
        isospec = chem_formula(elem).spectrum()
        # for each isotope plot line at mass:
        for mz in isospec:
            print(mz)
            masses_air_constituents.append(isospec[mz].mass)
            formulas_air_constituents.append(elem)




    
    path_combinatorial = r"C:\Users\kaho\Desktop\data\data_Empa\combinatorial_data"
    formulas_hydrocarbons_gs, formulas_hydrocarbons_rad, formulas_hydrocarbons, masses_hydrocarbons = get_chemical_possible_combinations_for_elem_list(elements = np.array(['C','H']), target_masses = np.arange(30, 300), radicals_allowed = True, path_to_files=path_combinatorial, recalculate=False)    
    
    print(masses_hydrocarbons.astype(int))

    # get exact masses of formulas

    calibrant_fragment_formulas_gs, calibrant_fragment_formulas_rad, calibrant_fragment_formulas, calibrant_fragment_masses = get_chemical_possible_combinations_for_elem_list(elements = np.array(['C','F']), target_masses = np.arange(30, 500), radicals_allowed=True,  path_to_files=path_combinatorial, recalculate=False)
    # get exact masses of formulas
    #masses_calibrant_fragments = [formula.isotope.mass for formula in calibrant_fragment_formulas]


    # plot the mass spectrum as before but with masses instead of indices as x-axis
    fig, ax = plt.subplots(1,1)
    # plot sum of all mass spectra
    #mass_sum_all = mass_sum.sum(axis=0)
    mass_sum_all = get_sum(data, sum_type= "mass", range_starts=left_rising_edge_index, range_ends=right_falling_edge_index)
    # shape of mass_sum_all: (number of ranges, number of mass bins)
    mass_sum_all_indici = np.arange(mass_sum_all.shape[0])

    plt.plot(mass_index_to_mass(mass_sum_all_indici)[0], mass_sum_all.sum(axis=1))
    # get indici of mass_sum_all
    #plt.plot(mz(mass_sum_all_indici), mass_sum_all)
    # plot the mass peaks
    plt.plot(index_to_mass(peaks[0], p1=p1, p2=p2), peaks[1]['prominences'], "x")
    # plot the masses of spikes
    plt.plot(masses, np.zeros(len(masses)), "x")
    # plot vertical lines at the masses of the formulas
    #plt.plot(masses_res2, np.zeros(len(masses_res2)), "x")
    for elem in masses_hydrocarbons:
        plt.axvline(x = elem,  color="red", linestyle="--", linewidth=0.5)
    for elem in calibrant_fragment_masses:
        plt.axvline(x = elem,  color="green", linestyle="--", linewidth=0.5)
    for elem in masses_air_constituents:
        plt.axvline(x = elem,  color="blue", linestyle="--", linewidth=0.5)


    # def linear_interpolate_intensity(x_eval_1D: np.array, x_1D, y_1D):
    #     """
    #     Linear interpolate the intensity of the mass spectra at the masses of the formulas.
    #     Parameters:
    #     masses_eval_1D: 1D array of masses
    #     mass_sum_all_indici_1D: 1D array of indices of the mass spectra
    #     mass_sum_all_multiD: multi-dimensional array of mass spectra
    #     """
    #     if isinstance(mass_sum_all_multiD, np.ndarray) and mass_sum_all_multiD.ndim == 1:
    #         # If mass_sum_all_multiD is 1D, convert it to a 2D array with one row
    #         mass_sum_all_multiD = mass_sum_all_multiD.reshape(1, -1)

    #     # get the indices of the mass spectra that are the next lower and next higher than the masses
    #     lower_indices = np.searchsorted(index_to_mass(mass_sum_all_indici), masses_eval_1D) - 1
    #     upper_indices = np.searchsorted(index_to_mass(mass_sum_all_indici), masses_eval_1D)
    #     # get the intensities at the lower and higher indices
    #     mass_sum_all_index_lower = mass_sum_all_multiD[:, lower_indices]
    #     mass_sum_all_index_upper = mass_sum_all_multiD[:, upper_indices]
    #     # get the delta mass between the lower and higher indices
    #     delta_mass = index_to_mass(mass_sum_all_indici[upper_indices]) - index_to_mass(mass_sum_all_indici[lower_indices])
    #     # set delta_mass to 1 if it is zero to avoid division by zero
    #     delta_mass[delta_mass == 0] = 1
    #     # calculate the intensity by linear interpolation
    #     intensity_linear_interpolated = mass_sum_all_index_lower + (masses_eval_1D - index_to_mass(mass_sum_all_indici[lower_indices])) / delta_mass * (mass_sum_all_index_upper - mass_sum_all_index_lower)
    #     return intensity_linear_interpolated





    masses_for_alpinac_hc = []
    formulas_for_alpinac_hc = []
    intensities_for_alpinac_hc = []
    masses_for_alpinac_cal = []
    formulas_for_alpinac_cal = []
    intensities_for_alpinac_cal = []
    masses_for_alpinac_air = []
    formulas_for_alpinac_air = []
    intensities_for_alpinac_air = []

    range_i=0
    intensity_hcs, masses_hcs, res = calculate_interpolated_intensities(masses = masses_hydrocarbons,
                            formulas = formulas_hydrocarbons, 
                            mass_indici = mass_sum_all_indici, 
                            intensities = mass_sum_all.sum(axis=1),
                            params_initial = params_initial,
                            threshold = 0.0,
                            plot=True
                            )
    
    # get parameters from res
    p1_fit = res.x[0]
    p2_fit = res.x[1]


    calculate_interpolated_intensities(calibrant_fragment_masses, calibrant_fragment_formulas, mass_sum_all_indici,  mass_sum_all.sum(axis=1), 0.0)
    #calculate_interpolated_intensities(masses_air_constituents, formulas_air_constituents,  mass_sum_all_indici, mass_sum_all[range_i], 0)

    # join all intensities and masses
    masses_all = np.concatenate((masses_hydrocarbons, calibrant_fragment_masses, masses_air_constituents))
    formulas_all = np.concatenate((formulas_hydrocarbons, calibrant_fragment_formulas, formulas_air_constituents))



    intensity_hcs, masses_hcs, res = calculate_interpolated_intensities(masses = masses_all,
                            formulas = formulas_all, 
                            mass_indici = mass_sum_all_indici, 
                            intensities = mass_sum_all.sum(axis=1),
                            params_initial = params_initial,
                            threshold = 0.0,
                            index_to_mass = index_to_mass,
                            plot=True


                            )

    print(res.fun)
    fit_par = res.x



    # CALIBRATION SPIKES

    # plot masses_all as vlines
    # Plot results
        # plot calibrant peaks in new calibration
    left_rising_edge_index, right_falling_edge_index, _, _, _, _ = get_calibrant_peaks(data_EI, start_or_end="start", plot = True)
    mass_sum_start = get_sum(data_EI, sum_type= "mass", range_starts=left_rising_edge_index, range_ends=right_falling_edge_index)
    left_rising_edge_index_end, right_falling_edge_index_end, _, _, _, _ = get_calibrant_peaks(data_EI, start_or_end="end", plot = True)
    mass_sum_end = get_sum(data_EI, sum_type= "mass", range_starts=left_rising_edge_index_end, range_ends=right_falling_edge_index_end)    
    formulas_spikes, masses_spikes = get_mass(spikes)

    # calibrate with calibrant peaks only:
    intensity_cal, masses_cal, res_cal = calculate_interpolated_intensities(masses = masses_all,
                            formulas = calibrant_fragment_formulas,
                            mass_indici = mass_sum_all_indici,
                            intensities = mass_sum_start.sum(axis=1),
                            params_initial = params_initial,
                            threshold = 0.0,
                            index_to_mass = index_to_mass,
                            plot=True
                            )
    fit_par_cal_start = res_cal.x

    # calibrate with calibrant peaks only:
    intensity_cal, masses_cal, res_cal = calculate_interpolated_intensities(masses = masses_all,
                            formulas = calibrant_fragment_formulas,
                            mass_indici = mass_sum_all_indici,
                            intensities = mass_sum_end.sum(axis=1),
                            params_initial = params_initial,
                            threshold = 0.0,
                            index_to_mass = index_to_mass,
                            plot=True
                            )
    fit_par_cal_end = res_cal.x


    #fit_par = fit_par_cal


    fig, ax = plt.subplots(1, 1)
    plt.plot(index_to_mass(mass_sum_all_indici, p1=fit_par[0], p2=fit_par[1]), mass_sum_start.sum(axis=1)) 
    plt.plot(index_to_mass(mass_sum_all_indici, p1=fit_par[0], p2=fit_par[1]), mass_sum_end.sum(axis=1))
    # add vertical lines at the masses of the calibrant peaks
    for i, mass in enumerate(masses_spikes):
        plt.axvline(x = mass,  color="red", linestyle="--", linewidth=0.5  )
        # label each line with formulas_spikes
        plt.text(mass, 0.1, formulas_spikes[i], rotation=90, fontsize=6)
    # add vertical lines at the masses of the formulas of calibrant_fragment_masses
    for i, mass in enumerate(calibrant_fragment_masses):
        plt.axvline(x = mass,  color="green", linestyle="--", linewidth=0.5  )
        # label each line with formulas_spikes
        plt.text(mass, 0.1, calibrant_fragment_formulas[i], rotation=90, fontsize=6)    


    # CI_data

    mass_sum_seg0 = get_sum(data = data_CI, sum_type= "mass", range_starts=np.array(0), range_ends=np.array(10000)   ) 
    mass_sum_seg_ind = np.arange(mass_sum_seg0.shape[0])    
    formulas_spikes, masses_spikes = get_mass(spikes)

    # calibrate with calibrant peaks only:
    intensity_cal, masses_cal, res_cal = calculate_interpolated_intensities(masses = masses_all,
                            formulas = calibrant_fragment_formulas,
                            mass_indici = mass_sum_all_indici,
                            intensities = mass_sum_start.sum(axis=1),
                            params_initial = params_initial,
                            threshold = 0.0,
                            index_to_mass = index_to_mass,
                            plot=True
                            )
    fit_par_cal_CI_seg = res_cal.x



    fig, ax = plt.subplots(1, 1)
    plt.plot(index_to_mass(mass_sum_seg_ind, p1=fit_par_cal_CI_seg[0], p2=fit_par_cal_CI_seg[1]), mass_sum_seg0.sum(axis=1), color = "black") 
    plt.plot(index_to_mass(mass_sum_all_indici, p1=fit_par[0], p2=fit_par[1]), mass_sum_end.sum(axis=1), color = "orange")
    # add vertical lines at the masses of the calibrant peaks
    for i, mass in enumerate(masses_spikes):
        plt.axvline(x = mass,  color="red", linestyle="--", linewidth=0.5  )
        # label each line with formulas_spikes
        plt.text(mass, 0.1, formulas_spikes[i], rotation=90, fontsize=6)
    # add vertical lines at the masses of the formulas of calibrant_fragment_masses
    for i, mass in enumerate(calibrant_fragment_masses):
        plt.axvline(x = mass,  color="green", linestyle="--", linewidth=0.5  )
        # label each line with formulas_spikes
        plt.text(mass, 0.1, calibrant_fragment_formulas[i], rotation=90, fontsize=6)  

    #TODO: get ranges of highest RT peaks, exclude them from data_CI and then do the calibration again
        # ranges with highest RT peaks
        #plot CI RT plot
        fig, ax = plt.subplots(1,1)
        plt.plot(data_CI.sum(axis=1))

    range_starts, range_ends = get_ranges_for_summing(data_CI, sum_type="RT", threshold_median_multiples=4, thresh_min_range_len= 3)
        # get the indices of the ranges
    range_starts_ind = 0

        #plt.plot(masses_hcs, intensity_hcs, "x")


# write an optimization routin which adapts the parameters of the index_to_mass function to obtain highest number of peaks over threshold and highest sum of intensities over threshold



    # # sum intensities of data at the masses of the formulas
    # masses_for_alpinac_hc = []
    # formulas_for_alpinac_hc = []
    # intensities_for_alpinac_hc = []
    # masses_for_alpinac_cal = []
    # formulas_for_alpinac_cal = []
    # intensities_for_alpinac_cal = []
    # masses_for_alpinac_air = []
    # formulas_for_alpinac_air = []
    # intensities_for_alpinac_air = []
    # for i, mass in enumerate(masses_hydrocarbons):
    #     # get the index of the mass in the mass_sum_all_indici array
    #     mass_index = np.where(mz(mass_sum_all_indici) == mass)
    #     # if found, calculate intensity:
    #     if not mass_index:
    #         # get the intensity at that index
    #         intensity_linear_interpolated = mass_sum_all[:, mass_index[0][0]]
    #     else:
    #         # if no index directly matched do linear interpolation of next two indices to get the intensity:
    #         # get the index of the next lower mass:
    #         lower_index = np.where(mz(mass_sum_all_indici) < mass)[0][-1]
    #         # get the index of the next higher mass:
    #         upper_index = np.where(mz(mass_sum_all_indici) > mass)[0][0]
    #         # get the intensity at the lower index
    #         mass_sum_all_index_lower = mass_sum_all[:, lower_index]
    #         # get the intensity at the higher index
    #         mass_sum_all_index_upper = mass_sum_all[:, upper_index]
    #         # calculate the intensity by linear interpolation
    #         intensity_linear_interpolated = mass_sum_all_index_lower + (mass - mz(mass_sum_all_indici[lower_index]))/(mz(mass_sum_all_indici[upper_index]) - mz(mass_sum_all_indici[lower_index])) * (mass_sum_all_index_upper - mass_sum_all_index_lower)
    #     # sum the intensities
    #     intensity_linear_interpolated_sum = intensity_linear_interpolated.sum()
    #     # if the sum is larger than 0.5, append to masses_res2
    #     if intensity_linear_interpolated_sum > 0.05:
    #         masses_for_alpinac_hc.append(mass)
    #         formulas_for_alpinac_hc.append(formulas_hydrocarbons[i])
    #         intensities_for_alpinac_hc.append(intensity_linear_interpolated_sum)
    # # same for calibrant fragments
    # for i, mass in enumerate(masses_calibrant_fragments):
    #     mass_index = np.where(mz(mass_sum_all_indici) == mass)
    #     # if found, calculate intensity:
    #     if not mass_index:
    #         intensity_linear_interpolated = mass_sum_all[:, mass_index[0][0]]
    #     else:
    #         # get the intensity at that index
    #         lower_index = np.where(mz(mass_sum_all_indici) < mass)[0][-1]
    #         # get the index of the next higher mass:
    #         upper_index = np.where(mz(mass_sum_all_indici) > mass)[0][0]
    #         # get the intensity at the lower index
    #         mass_sum_all_index_lower = mass_sum_all[:, lower_index]
    #         # get the intensity at the higher index
    #         mass_sum_all_index_upper = mass_sum_all[:, upper_index]
    #         # calculate the intensity by linear interpolation
    #         intensity_linear_interpolated = mass_sum_all_index_lower + (mass - mz(mass_sum_all_indici[lower_index]))/(mz(mass_sum_all_indici[upper_index]) - mz(mass_sum_all_indici[lower_index])) * (mass_sum_all_index_upper - mass_sum_all_index_lower)
    #     # sum the intensities
    #     intensity_linear_interpolated_sum = intensity_linear_interpolated.sum()
    #     # if the sum is larger than 0.5, append to masses_res2
    #     if intensity_linear_interpolated_sum > 0.05:
    #         masses_for_alpinac_cal.append(mass)
    #         formulas_for_alpinac_cal.append(calibrant_fragment_formulas[i])
    #         intensities_for_alpinac_cal.append(intensity_linear_interpolated_sum)
     
    # # same for air constituents
    # for i, mass in enumerate(masses_air_constitutents):
    #     mass_index = np.where(mz(mass_sum_all_indici) == mass)
    #     # if found, calculate intensity:
    #     if not mass_index:
    #         intensity_linear_interpolated = mass_sum_all[:, mass_index[0][0]]
    #     else:
    #         lower_index = np.where(mz(mass_sum_all_indici) < mass)[0][-1]
    #         # get the index of the next higher mass:
    #         upper_index = np.where(mz(mass_sum_all_indici) > mass)[0][0]
    #         # get the intensity at the lower index
    #         mass_sum_all_index_lower = mass_sum_all[:, lower_index]
    #         # get the intensity at the higher index
    #         mass_sum_all_index_upper = mass_sum_all[:, upper_index]
    #         # calculate the intensity by linear interpolation
    #         intensity_linear_interpolated = mass_sum_all_index_lower + (mass - mz(mass_sum_all_indici[lower_index]))/(mz(mass_sum_all_indici[upper_index]) - mz(mass_sum_all_indici[lower_index])) * (mass_sum_all_index_upper - mass_sum_all_index_lower)
    #     # sum the intensities
    #     intensity_linear_interpolated_sum = intensity_linear_interpolated.sum()
    #     # if the sum is larger than 0.5, append to masses_res2
    #     if intensity_linear_interpolated_sum > 0.05:
    #         masses_for_alpinac_air.append(mass)
    #         formulas_for_alpinac_air.append(air_constitutents[i])
    #         intensities_for_alpinac_air.append(intensity_linear_interpolated_sum)


    # sort the masses and formulas by mass
    masses_for_alpinac = np.array(masses_for_alpinac)
    formulas_for_alpinac = np.array(formulas_for_alpinac)
    intensities_for_alpinac = np.array(intensities_for_alpinac)

    # combine all arrays:
    masses_for_alpinac = np.concatenate((masses_for_alpinac_hc, masses_for_alpinac_cal, masses_for_alpinac_air))
    formulas_for_alpinac = np.concatenate((formulas_for_alpinac_hc, formulas_for_alpinac_cal, formulas_for_alpinac_air))
    intensities_for_alpinac = np.concatenate((intensities_for_alpinac_hc, intensities_for_alpinac_cal, intensities_for_alpinac_air))







    # # 
    # # Give a list of all files ending with .txt in the directory:
    # import os
    # dir = r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303"
    # files_h5 = [f for f in os.listdir(dir) if f.endswith(".h5")]
    # files_mc = [f for f in os.listdir(dir) if f.endswith(".mc.txt")]
    # # get all files that do not have a corresponding .mc.txt file
    # files_h5_no_mc = []
    # for file_h5 in files_h5:
    #     if not file_h5[:-2] + "mc.txt" in files_mc:
    #         files_h5_no_mc.append(file_h5)

    # len(files_h5_no_mc)
    # len(files_mc)
    # len(files_h5)



    # # read all files and append to a list
    # # copy all files of files_h5_no_mc to a new directory:
    # import shutil
    # dir_new = r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\files_h5_no_mc"
    # # create new directory
    # #os.mkdir(dir_new)
    # for file_h5 in files_h5_no_mc:
    #     shutil.copy(os.path.join(dir, file_h5[:-2] + "mc.txt"), dir_new)

    # # copy all files where it worked to a new directory:
    # dir_new = r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\files_h5_mc"
    # # create new directory
    # #os.mkdir(dir_new)
    # for file_h5 in files_h5:
    #     if file_h5[:-2] + "mc.txt" in files_mc:
    #         # save in new directory
    #         shutil.copy(os.path.join(dir, file_h5[:-2] + "mc.txt"), dir_new)


    
    # 




    print(masses_res2)

    print(formulas_res2)

    # add to plot
    #plt.plot(peaks, mass_sum_all[peaks], "x")




    # select range for first fit, 0 to second flank
    fit_start = 5
    fit_end = flank_indici[1]
    # fit a sigmoidal exponential to the data
    xdata = rt_sum_threshold_index_interval_index[fit_start:fit_end]
    ydata = rt_sum_threshold_index_interval[fit_start:fit_end]
    fig, ax = plt.subplots(1,1)
    plt.plot(xdata, ydata, color="blue", linewidth=0.5, label="data")

    # find value where the data is 50% of the maximum
    half_max_index = np.where(ydata < half_max)[0][0]
    half_max_value = ydata_rescaled[half_max_index]



    # # plot a sigmoidal exponential
    # def approx_heavy_expo_decay(t, alpha, beta, tpeak, A):
    #     return A/(1+np.exp(-alpha*(t-tpeak)))* np.exp(-beta * (t-tpeak))

    # # Fit the data to the function
    # popt, pcov = curve_fit(f = approx_heavy_expo_decay, xdata = xdata, ydata = ydata_rescaled, p0 = [20000, 0.5, 610, 1])

    # # Print the fitted parameters
    # print("A = {:.3f}".format(popt[0]))
    # print("a = {:.3f}".format(popt[1]))
    # print("b = {:.3f}".format(popt[2]))
    # print("xstep = {:.3f}".format(popt[3]))
    # #popt = [2000, 0.5, 610, 1]
    # plt.plot(xdata, approx_heavy_expo_decay(xdata, *popt), color="red", linestyle="-", linewidth=0.5, label="fit")