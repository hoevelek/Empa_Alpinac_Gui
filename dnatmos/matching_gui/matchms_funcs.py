# -*- coding: utf-8 -*-
"""
Created on Fri May 20 19:08:22 2022

@author: kaho
"""
from __future__ import annotations

# todo convert Myriam file to mgf
# what is inchikey
# what is derive inchi from smiles

# system
from pathlib import Path
import os
import sys
import logging

# data/math
import numpy as np
import pandas as pd
import re #regex
from itertools import chain
import csv

# plotting
from matplotlib import pyplot as plt, rcParams
import matplotlib.patches as mpl_patches #for getting label handles
import matplotlib.ticker as mticker
from matplotlib import axes
import PIL.Image
import PIL.IcoImagePlugin
import PIL#plotting pngs


# matchms packages
import matchms as ms
from matchms.similarity import ModifiedCosine
from rdkit.Chem import Draw
from rdkit import Chem
from matchms.similarity import CosineGreedy
from matchms import calculate_scores
from matchms.filtering import select_by_mz
from matchms.filtering import select_by_intensity
from matchms.filtering import normalize_intensities
from matchms.filtering import add_precursor_mz
from matchms.filtering import harmonize_undefined_smiles
from matchms.filtering import harmonize_undefined_inchikey
from matchms.filtering import harmonize_undefined_inchi
from matchms.filtering import derive_inchi_from_smiles
from matchms.filtering import derive_smiles_from_inchi
from matchms.filtering import derive_inchikey_from_inchi
from matchms.filtering import repair_inchi_inchikey_smiles
from matchms.filtering import default_filters
# to load from a specific file format
from matchms.importing import load_from_mgf
from matchms.exporting import save_as_json
from matchms.exporting import save_as_mgf
from matchms.plotting.spectrum_plots import plot_spectrum

# for StickSpectra
from lmfit import Minimizer
from lmfit.parameter import Parameters


# files
import pickle

# scientific
from periodictable import elements
# for calculating molecular formula from smiles
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
# to convert to latex formula #x='C'
from pyvalem.formula import Formula
# for calculation of isotopic pure mass and isotopic spectra
from molmass import Formula as chem_formula



from typing import Optional

# logger
from logging import getLogger, StreamHandler, Formatter, FileHandler
logger = getLogger(__name__)
logger.setLevel(logging.DEBUG)


def metadata_processing(spectrum):
    spectrum = default_filters(spectrum)
    spectrum = repair_inchi_inchikey_smiles(spectrum)
    spectrum = derive_inchi_from_smiles(spectrum)
    spectrum = derive_smiles_from_inchi(spectrum)
    spectrum = derive_inchikey_from_inchi(spectrum)
    spectrum = harmonize_undefined_smiles(spectrum)
    spectrum = harmonize_undefined_inchi(spectrum)
    spectrum = harmonize_undefined_inchikey(spectrum)
    spectrum = add_precursor_mz(spectrum)
    return spectrum

def peak_processing(spectrum):
    spectrum = normalize_intensities(spectrum)
    spectrum = select_by_intensity(spectrum, intensity_from=0.01)
    spectrum = select_by_mz(spectrum, mz_from=10, mz_to=1000)
    return spectrum

def print_extra_info_in_label(spectrum: ms.Spectrum, ax, no_plt:int = 0):
    """prints extra (header) information of spectra in a label"""
    # work around to place text at best position, needed if loc = 'best'
    if not 'smiles' in spectrum.metadata:
        logging.warn(r"No smiles code for {}.replace(spectrum)")
        leg = None
        return(leg)
    smiles_code = spectrum.get("smiles")
    handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white",
                                     lw=0, alpha=0)] * 4
    mol = Chem.MolFromSmiles(smiles_code)
    formula = CalcMolFormula(mol)
    hetatms = {atom.GetSymbol() for atom in mol.GetAtoms()} #TODO
    latex_formula = formula
    latex_formula = Formula(formula).latex
    
    labels = []
    # labels.append(r'smiles: {}'.format(spectrum.get("smiles")))
    # labels.append(r'name: {}'.format(spectrum.get("name")))
    labels.append(r'formula: ${}$'.format(latex_formula))
    labels.append(r'mass/u: {}'.format(spectrum.get("precursor_mz")))
    # create the legend, supressing the blank space of the empty line symbol and the
    # padding between symbol and label by setting handlelenght and handletextpad
    if no_plt > 0:
        title_i = r'Proposed molecule #{}'.format(no_plt)
    else:
        title_i = 'Requested molecule'
    if no_plt == 0:
        leg = ax.legend(handles, labels, loc='upper right', fontsize='x-small',
                        fancybox=True, framealpha=0.7,
                        handlelength=0, handletextpad=0, title=title_i)
    else:
        leg = ax.legend(handles, labels, loc='lower right', fontsize='x-small',
                           fancybox=True, framealpha=0.7,
                           handlelength=0, handletextpad=0, title=title_i)
    plt.setp(leg.get_title(), fontsize='x-small',color=r'C{}'.format(no_plt))
    plt.setp(leg.get_texts(), color=r'C{}'.format(no_plt))
    # add to make sure it is not overwritten vom second legend
    plt.gca().add_artist(leg)
    return(leg)

def score_matching(query_spectrum: ms.Spectrum, 
                   ref_spectra: ms.Spectrum | list[ms.Spectrum] = [],
                   tolerance_level:float = 0.005, min_match_peaks: int = 5,
                   no_best_matches: int = 10) -> tuple[list[ms.Spectrum], list]:
    """ calculates the matching (number of commun peaks and cosine greedy
    scores within a tolerance mass shift (shift in m/z)) and returns an
    object with first {no_best_matches} best matches, a list
    with a (matchms spectrum, cosine score, and #common peaks)"""
    similarity_measure = CosineGreedy(tolerance=tolerance_level)
    # Asymmetric, compare a certain spectrum to a library (no mass shift considered)
    if isinstance(query_spectrum, ms.Spectrum):
        query_spectrum = [query_spectrum] #convert to list, as req. for calculate scores
    else:
        logger.info ("No valid query spectrum input. Note: a list of spectra is not permitted for 1D search.")
    if isinstance(ref_spectra, ms.Spectrum):
        ref_spectra = [ref_spectra]
    
    scores1D = calculate_scores(references=ref_spectra, queries=query_spectrum,
                                similarity_function=similarity_measure, is_symmetric=False)
    # scores = scores1D.scores['score']
    # no_matching_peaks = scores1D.scores['matches']

    
    
    #for (reference, query, score) in scores1D:
    #    logger.info(f"Cosine score between {reference.get('id')} and {query.get('id')}" +
    #          f" is {score['score']:.2f} with {score['matches']} matched peaks")
    # number of matching peaks
    # scores1D.scores[:len(ref_spectra),:len(ref_spectra)]["matches"]
    # cosine score
    # scores1D.scores[:len(ref_spectra), :len(ref_spectra)]["score"]
    # require a min number of {no_best_matches} fitting peaks
    # TODO check if not ref_sepctra
    sorted_matches = scores1D.scores_by_query(query_spectrum, sort=True)
    #logger.info([x[1]["score"].round(3) for x in sorted_matches])
    best_matches = [x for x in sorted_matches if x[1]
                    ["matches"] >= min_match_peaks][:no_best_matches]

    #best_matches = [*query_spectrum,*query_spectrum]
    scores = [0 for i in best_matches]   
    best_matches_spectra = [el[0] for el in best_matches]
    best_matches_score = [el[1]["score"] for el in best_matches]
    best_matches_no_matching_peaks = [el[1]["matches"] for el in best_matches]
    logger.info(best_matches_spectra)
    logger.info(best_matches_score)
    return best_matches_spectra, best_matches_score, best_matches_no_matching_peaks

def plot_mass_spectra(ax, spectrum: ms.Spectrum, 
                      spectrum2: ms.Spectrum | list[ms.Spectrum] = [],
                      no_plt_legend: int = 1, **kwargs):
    """Plots a first, e.g. the experimental unknown spectra against the a list
    of other spectra to compare (upside down)
    Arguments:
        ax: matplotlib axis
        spectrum: matchms.Spectrum
        spectrum2: matchms.Spectrum or list of matchms.Spectrum
        no_plt_legend: int, if 0, no legend is plotted
        **kwargs: additional arguments passed to matplotlib.pyplot.plot
    """

    # get which peak label is shown
    if 'peak_labels' in kwargs: 
        labels_peak = kwargs['peak_labels']                 
    else:
        labels_peak = 'formulas'

    # get figsize if given
    if 'figsize' in kwargs:
        figsize = kwargs['figsize']
    else:
        fig_width, fig_height = ax.get_figure().get_size_inches()
        figsize = (fig_width, fig_height)

    if 'col' in kwargs:
        col = kwargs['col']
    else:
        col = 'C0'

    # check if figsize is correct
    if not isinstance(figsize, tuple) or len(figsize) != 2 or \
    not all(isinstance(v, (int, float)) and v > 0 for v in figsize):
        # If figsize is invalid, set it to a default value
        figsize = (10, 8)
        print("Invalid figsize argument. Setting to default value (10, 8).")

    # get label if given
    if 'label' in kwargs:
        label_stem = kwargs['label']
        # if label is a list, check if it has the correct length
        if isinstance(label_stem, list):
            if len(label_stem) != 2:
                raise ValueError("The input argument 'label' must have a length of 2. First is for the stmlines, second for verticval lines")
    else:
        label_stem = ["theor.", "meas."]
    if isinstance(spectrum, AlpinacData):
        spectrum.plot_stem(ax, **kwargs)
    elif not len(spectrum.peaks.mz):
        # Check empty spectrum
        pass
    else:
        markerline, stemline, baseline = ax.stem(
            spectrum.peaks.mz, spectrum.peaks.intensities)
        plt.setp(stemline, linewidth=1.0,color=col)
        plt.setp(baseline, linewidth=0.5, color="black")
        plt.setp(markerline, markersize=1.5,color=col)
        ax.set_xlabel(r"$m\//z$ (Da)")
        ax.set_ylabel(r"Norm. Intensity (arb. u.)")
    # ax.tick_params(axis='y', colors='C0')
    # convert ms spectrum to a list if not in a list
    if isinstance(spectrum2, ms.Spectrum):
        spectrum2 = [spectrum2]
    # print legend for queried spectrum
    print_extra_info_in_label(spectrum, ax, 0)
    # add to make sure it is not overwritten vom second legend
    # plt.gca().add_artist(leg1)
    labels = ['Queried spectrum']

        # print second info 
    if spectrum2:
        ii = 0
        for spec2 in spectrum2: #spec2=spectrum2[0]
            if ii==0:
                ax.set(ylim=(-1.2, 1.2))
            ii += 1
            if isinstance(spec2, AlpinacData):
                spectrum.plot_stem(ax, **kwargs)
            else:
                markerline, stemline, baseline = ax.stem(
                spec2.peaks.mz, -spec2.peaks.intensities, label=r'#{}'.format(ii))
                plt.setp(baseline, linewidth=0.5, color="black")
                plt.setp(markerline, markersize=1.5, color=r'C{}'.format(ii))
                plt.setp(stemline, linewidth=1.0, color=r'C{}'.format(ii))
            labels.append(r'Mol. #: {}'.format(ii))
        print_extra_info_in_label(spectrum2[int(no_plt_legend-1)], ax, no_plt_legend)
    # add second axis with same scaling
    ax.yaxis.set_major_formatter(lambda x, pos: f'{abs(x):g}')

    
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),
          ncol=5, fancybox=True, shadow=True, fontsize = 'x-small')
    
def get_header_info(
        # return header info of a MassBank-mfg file for plotting 
        spectra:list[ms.Spectrum],
        metadata:str | list[str]=['pepmass', 'charge', 'mslevel', 
                                        'source_instrument', 'seq', 'ionmode',
                                        'organism', 'smiles', 'inchi',
                                        'inchiaux', 'pubmed', 'libraryquality',
                                        'scans', 'file_name', 'compound_name',
                                        'principal_investigator', 'data_collector',
                                        'submit_user', 'spectrum_id', 'precursor_mz',
                                        'adduct', 'inchikey'],
        print_to_table:bool = True,
        print_mol_struct:bool = True
    ) -> pd.DataFrame:
    """return header info of a MassBank-mfg file for plotting """
    """Return the metadata of the given spectrum for the given metadata parameters."""
    if isinstance(spectra, ms.Spectrum):
        spectra = [spectra]
    if isinstance(metadata, str):
        metadata = [metadata]

    dic = {
        spectrum: [
            spectrum.metadata.get(key, None)
            for key in metadata
        ]
        for spectrum in spectra
    }
    
    df =  pd.DataFrame.from_dict(dic, orient="index", columns=metadata)
    return df


class AlpinacData(ms.Spectrum):
    # TODO: save into loadable format (mgf), update .plot, mean spectra, std spectra, Name to Formula
    # Missing .uncertainties, .formulas, .dev_from_theo ..., .name, ..
    # transfer keynames to those used by matchms
    # alpinacData.input.spectra
    # alpinacData.output.spectra
    # test if works for input data

    uncertainties: np.array #initializes uncertainities as global variable
    def __init__(
        self,
        mz: np.array,
        intensities: np.array,
        metadata: Optional[dict] = None,
        metadata_harmonization: bool = True
    ):
        super().__init__(mz, intensities, metadata, metadata_harmonization) #super: inheritance from parent class
        # Store the uncertainties
        self._uncertainties = None #only internal
    

    #Is this needed?
    #@property #property: called as x.property instead of x.property()
    @classmethod 
    def _get_property_from_metadata(cls, _property_name) -> np.ndarray:
        """ get properties from metadata attribute if not already setted."""
        if getattr(cls, _property_name) is None: #if internal property uncertainty not set yet, it will set it
            if _property_name in cls.metadata:
                # Read from the metadata
                metadata_property = cls.metadata[_property_name]
                # TODO: convert the metadata str
                setattr(cls, name = _property_name, value = metadata_property)
            raise ValueError(f"{_property_name} not set in {cls}.")
        return getattr(cls, _property_name)

    #@_set_property_to_metadata.setter # sets properties, so that x.uncertainties = y works
    @classmethod 
    def _set_property_to_metadata(cls, _property_name, _value: np.array):
        """Set the value of the properties as string converted numpy array in the metadata."""
        if not isinstance(_value, np.ndarray):
            raise TypeError(f"{_property_name} must be a np.array")
        cls._property_name = _value
        # Also set in the metadata
        setattr(cls, _property_name, _value)
        #cls.set("delta_mz_ppm_meta", _value) #TODO does not work
        print("set metadata")
        #AlpinacData.set[_property_name] = str(list(_value))
        #return _value

    #@classmethod #decorator    
    #TODO add all properties which must be handed over to metadata for database search 
    @property #property: called as x.property instead of x.property()
    def delta_mz_ppm_meta(self) -> np.ndarray:
        self._get_property_from_metadata("delta_mz_ppm_meta")
   
    @delta_mz_ppm_meta.setter #should change in metadata and in attribute
    def delta_mz_ppm_meta(self, value):
        val = self._set_property_to_metadata("delta_mz_ppm_meta", value)
        self.set("delta_mz_ppm_meta", value)

    @classmethod
    def _clean_str_list(cls, str_list, repl_dict): #general function TODO -> in module tools
    # remove unsuited characters to make, e.g. a conversion of a column name to a variable possible
        if isinstance(str_list, str):
            for key in repl_dict.keys(): #key = ','
                return_list = str_list.replace(key, repl_dict[key])
        if isinstance(str_list, list):
            return_list = []
            for elem in str_list:
                for key in repl_dict.keys(): #key = ','
                    elem = elem.replace(key, repl_dict[key])
                return_list.append(elem.strip())
        return return_list
    
    
    @classmethod #decorator  
    def _autodetect_output_header(cls, path): 
        """extracts scalar headers and header table. It assumes first data line the one with a numeric value and would not work if a headerline starts with a number."""
        header_dict={} #for info in the format key:value
        flag_header_data = 0 # to know if we are in an key or value line (iterating)
        repl_dict_vars = {':':'',' ':''} #strings to replace in variables (replace key by value)
        repl_dict_keys = {'.':'',',':'','(':'in ','[':'in ',' ':'_','%':'percent',')':'',']':'','?':'',';':''} # strings to replace in variable names/ column names (replace key by value)
        no_block_cols = 1
                
        with open(path, 'r') as f:
            reader = csv.reader(f, delimiter = '\t')
            for index, row in enumerate(reader):
                if row[0].startswith('#'): #skip seperator lines
                    continue
                if row[0].startswith('Largest'): # skip info line
                    continue

                # adapt row such that it can be seperated into words, adaption for new formate
                if len(row)==1:
                    row_content = row[0].strip().replace('Likelihood indicator (%)','Likelihood_indicator_in_percent;').replace('Molecular ion? ', 'Molecular ion?;').replace('DBE', 'DBE;').replace('Formula', 'est_precursor_formula;').replace('+:', '+')
                    #splitting header line into words
                    row = row_content.split(';')
                    # if it is still len one than it is maybe header info seperated by a :
                    # if one word header line of format "key: value" (used for header info, not seperated by \t compared to data). 
                    if len(row)==1: #equiv. to if line with not vectorized info # new format by Blendi: does not work anymore! All lines have random white spaces as seperators, as well as the variable names ...
                        rsp_h = row[0].split(':')
                        rsp = [x for x in rsp_h if x.strip()]
                        if len(rsp)==2:
                            # remove unsuitable chars and add to header name dictionary
                            headervar = AlpinacData._clean_str_list(rsp[0],repl_dict_keys)
                            header_dict[headervar] = AlpinacData._clean_str_list(rsp[1],repl_dict_vars)
                        else: #must be data
                            row = row[0].split()

                print(row)
                if row[0].replace('.','',1).isdigit(): # now the numeric data would follow
                    print("header has so many lines:", index-1)
                    break
                              
                # if haeder table format (multiple words): first occurent row is used for colnames of header table
                # does not work anymore for newer files, where everything is sepearted by a random number of spaces
                if len(row)>1:
                    # strip white spaces 
                    row = [x.strip() for x in row if x.strip()]
                    if flag_header_data==0: #initiate np.array with colnames
                        block = np.array(AlpinacData._clean_str_list(row, repl_dict_keys))
                        flag_header_data=1
                        no_block_cols = len(row)
                    elif len(row)==no_block_cols:
                        block = np.vstack([block, AlpinacData._clean_str_list(row, repl_dict_vars)])
                        flag_header_data+=1
                    else:
                        pass
                print(row)
                        #print("block ended")
        # if table-like metadata this should contain more than just the colnames / only colnames could also mean data colnames are missinterpreted as metadata colnames
        if(flag_header_data>1):             
            df_b = pd.DataFrame(block)
            df_b = df_b.rename(columns=df_b.iloc[0]).drop(df_b.index[0]).reset_index(drop=True)
            for key in df_b.columns:
                header_dict.update({key: df_b[key].to_string(index=False).replace('\n',', ')})
                if key == "Formula": 
                    header_dict["est_precursor_formula"] = header_dict['Formula'] #prevent clash with fragment formula
                    del header_dict['Formula']
        return [header_dict, index-1] 

    def to_alpinac(self, path: Path):
        """Save this to a file for alpinac input/output data."""
        ...
    
    def to_database_mgf(self, path:Path):
        """ converts to mgf library. """
        # TODO convert only EI data to mgf, select what shoulb be in header and match required parameters """
        # TODO modify header TODO, what should be identifier if name, formula are unknown?
        save_as_mgf(self, path)
    
    @classmethod #decorator
    def from_alpinac(cls, path: Path, **kwargs) -> AlpinacData: 
        """this function should be a method of class AlpinacData and convert an alpinac generated txt file into a class AlpinacData-Spectrum.
        It converts input and output data into an ms.spectrum, adds non mass and intensity values as attributes, set scalar header data automatically by assigning the char before ":" to the value thereafter
        and sets vectorized header corresponding to the output of the scoring algorith by alpinac (basically a data table before the data table with single lines) to header values and attributes. 
        **kwargs:
            mass_return: str, default = 'exact_mass_in_m_per_z'
                'exact_mass_in_m_per_z' -> returns exact mass in m/z
                alternative: 'meas_mass_in_m_per_z' -> returns measured mass in m/z
        """
        # optional argum
        # check if optional argument mass_return is given
        if 'mass_return' in kwargs:
            mass_return = kwargs['mass_return']
        else:
            mass_return = 'exact_mass_in_m_per_z'

        if 'ionization_mode' in kwargs:
            ionization_mode = kwargs['ionization_mode']
        else:
            ionization_mode = 'EI'
        # TODO: see description above
        # process header
        return_header_detect = AlpinacData._autodetect_output_header(path) #TODO WHY does it not work
        last_line_header = return_header_detect[1]
        header_as_dict = return_header_detect[0]
        # write header info to metadata
        print("load spectrum")

        # read only header line in:
        header_line = pd.read_csv(path, sep = ';',header=last_line_header, nrows=0)

        if len(header_line.columns)>1: #new formate
            column_names = [x.strip() for x in header_line.columns if x.strip()]
            df = pd.read_csv(path, sep='\s+', names=column_names, skiprows=last_line_header+1)
        else: #old format
            df = pd.read_csv(path, sep = '\t',header=last_line_header)

        # bring column names to variable format
        # TODO use new function
        repl_list = {'.':'',',':'','(':'in ','[':'in ',')':'',']':'',' ':'_','%':'percent','/':'_per_'} # to remove everything which is a character not useful for variables
        for name in df.columns:
            old_name = name
            for key in repl_list.keys(): #key = ','
                name = name.replace(key, repl_list[key]).lower()
            df.rename(columns = {old_name:name}, inplace = True)
        print(f"file columns are {df.columns}")

        # use only df with ionization mode of interest
        df = df[df['ionization_mode']==ionization_mode]
        # dependent of if a key is available define spectrum/properties   
        #convert to spectrum
        #try if ouput data (first statement, specific keys) or input data (other specific data). Extendable to other file formates by adding elifs
        if 'exact_mass_in_m_per_z' in df.columns and 'percent_of_total_intensity' in df.columns:
            df = df.sort_values(by = 'exact_mass_in_m_per_z') #matchms class can only convert to spectrum if ordered by increasing mz values
            mz = df[mass_return].to_numpy()
            intensities = df['percent_of_total_intensity'].to_numpy()/100
        elif 'mass' in df.columns and 'area' in df.columns: # this is the case if you use an alpinac input file for converting
            df = df.sort_values(by = 'mass')
            mz = df['mass'].to_numpy()
            intensities = df['area'].to_numpy()/max(df['area'])
            mass_return = 'mass'
        else:
            errmess = r"Column names {} were not detected in data file. ".format(''.join([str(item)+" | " for item in df.columns]))+\
                "Expecting 'Exact_mass_in_m_per_z' or"+\
                "'mass' and 'percent_of_Total_Intensity' or  'area'. "+ \
                "Check file or add elif block for new file format"
            raise KeyError(errmess)
        
        new_spectrum = AlpinacData(mz = mz, intensities = intensities) # = AlpinacData()

        print("set metadata")
        for key in header_as_dict.keys():
            new_spectrum.set(key, header_as_dict[key])

        print("vector-type metadata")
        # TODO extend by what you want
        if 'mass_u_ppm' in df.columns:
            new_spectrum._set_property_to_metadata("delta_mz_ppm_meta", df['mass_u_ppm'].to_numpy())
            new_spectrum.set("delta_mz_ppm_meta", str(list(df['mass_u_ppm'].to_numpy())))
    
        print('set non-matchms-type data as attributes')
        blacklist_attributes = ['mass', mass_return, 'area', 'percent_of_total_intensity']
        for colname in [elem for elem in df.columns if elem not in blacklist_attributes]:
            setattr(new_spectrum, colname, df[colname])
        return new_spectrum
    
    @classmethod #decorator
    def from_alpinac_unidentified(cls, path: Path) -> AlpinacData: 
        """this function should be a method of class AlpinacData and convert an alpinac generated txt file into a class AlpinacData-Spectrum.
        It returns all peaks which are not(!) identified and returns its measured mass and intensity.
        It converts input and output data into an ms.spectrum, adds non mass and intensity values as attributes, set scalar header data automatically by assigning the char before ":" to the value thereafter
        and sets vectorized header corresponding to the output of the scoring algorith by alpinac (basically a data table before the data table with single lines) to header values and attributes. 
        """
        # TODO: see description above
        # process header
        return_header_detect = AlpinacData._autodetect_output_header(path) #TODO WHY does it not work
        last_line_header = return_header_detect[1]
        header_as_dict = return_header_detect[0]
        # write header info to metadata
        print("load spectrum")

        # read only header line in:
        header_line = pd.read_csv(path, sep = ';',header=last_line_header, nrows=0)

        if len(header_line.columns)>1: #new formate
            column_names = [x.strip() for x in header_line.columns if x.strip()]
            df = pd.read_csv(path, sep='\s+', names=column_names, skiprows=last_line_header+1)
        else: #old format
            df = pd.read_csv(path, sep = '\t',header=last_line_header)
        # bring column names to variable format
        # TODO use new function
        repl_list = {'.':'',',':'','(':'in ','[':'in ',')':'',']':'',' ':'_','%':'percent','/':'_per_'} # to remove everything which is a character not useful for variables
        for name in df.columns:
            old_name = name
            for key in repl_list.keys(): #key = ','
                name = name.replace(key, repl_list[key]).lower()
            df.rename(columns = {old_name:name}, inplace = True)
        print(f"file columns are {df.columns}")
        # use only df for which exact mass is 0
        df = df[df['exact_mass_in_m_per_z'] == 0]
        # if no data is left, return None
        # TODO: check if this is the best way to handle this

        # dependent of if a key is available define spectrum/properties   
        #convert to spectrum
        #try if ouput data (first statement, specific keys) or input data (other specific data). Extendable to other file formates by adding elifs
        if 'meas_mass_in_m_per_z' in df.columns and 'percent_of_total_intensity' in df.columns:
            df = df.sort_values(by = 'meas_mass_in_m_per_z') #matchms class can only convert to spectrum if ordered by increasing mz values
            mz = df['meas_mass_in_m_per_z'].to_numpy()
            intensities = df['percent_of_total_intensity'].to_numpy()/100
        else:
            errmess = r"Column names {} were not detected in data file. ".format(''.join([str(item)+" | " for item in df.columns]))+\
                "Expecting 'Exact_mass_in_m_per_z' or"+\
                "'mass' and 'percent_of_Total_Intensity' or  'area'. "+ \
                "Check file or add elif block for new file format"
            raise KeyError(errmess)
        
        new_spectrum = AlpinacData(mz = mz, intensities = intensities) # = AlpinacData()

        print("set metadata")
        for key in header_as_dict.keys():
            new_spectrum.set(key, header_as_dict[key])

        print("vector-type metadata")
        # TODO extend by what you want
        if 'mass_u_ppm' in df.columns:
            new_spectrum._set_property_to_metadata("delta_mz_ppm_meta", df['mass_u_ppm'].to_numpy())
            new_spectrum.set("delta_mz_ppm_meta", str(list(df['mass_u_ppm'].to_numpy())))
    
        print('set non-matchms-type data as attributes')
        blacklist_attributes = ['mass','exact_mass_in_m_per_z','area', 'percent_of_total_intensity']
        for colname in [elem for elem in df.columns if elem not in blacklist_attributes]:
            setattr(new_spectrum, colname, df[colname])
        return new_spectrum
    
    # this function needed to be renamed from plot to plot_stem because of a clash of names with the matchms class
    def plot_stem(self, ax, **kwargs): #TODO df to alpinac class
        """ plot a stick spectra of the measured (input data) and ALPINAC-calculated data (output data). 
        In the later case it assigns peaks to fragments and plots measured data as a vertical line
        Arguments:
        `ax: Matplotlib axes object
        Optional arguments:
        - figsize: tuple of length 2, size of the figure
        - col: color of the stem plot
        - label: list of length 2, labels for the stem plot and the vertical lines
        - label_peaks: str, either "formulas" or "numeric" or None, if "formulas" the formula of the fragment is plotted, if "masses" the mass of the fragment is plotted, default is "formulas"
        """
        # check if ax is a matplotlib axes object
        if not isinstance(ax, plt.Axes):
            raise TypeError("The input argument 'ax' must be a Matplotlib axes object.")

        # check if optional argument colr (for plotting stick spectrum) is given
        if 'col' in kwargs:
            col_stem = kwargs['col']
        else:
            col_stem = 'C0'
        # get figsize if given
        if 'figsize' in kwargs:
            figsize = kwargs['figsize']
        else:
            fig_width, fig_height = ax.get_figure().get_size_inches()
            figsize = (fig_width, fig_height)

        if "label_peaks" in kwargs:
            label_peaks = kwargs["label_peaks"]
        else:
            label_peaks = "formulas"

        # check if figsize is correct
        if not isinstance(figsize, tuple) or len(figsize) != 2 or \
        not all(isinstance(v, (int, float)) and v > 0 for v in figsize):
            # If figsize is invalid, set it to a default value
            figsize = (10, 8)
            print("Invalid figsize argument. Setting to default value (10, 8).")

        # get label if given
        if 'label' in kwargs:
            label_stem = kwargs['label']
            # if label is a list, check if it has the correct length
            if isinstance(label_stem, list):
                if len(label_stem) != 2:
                    raise ValueError("The input argument 'label' must have a length of 2. First is for the stmlines, second for verticval lines")
        else:
            label_stem = ["theor.", "meas."]


        # get fragments name
        substance = "unknown"
        if hasattr(self, "formula") and label_peaks == "formulas":
                formulas = self.formula
        elif not hasattr(self, "formula") and label_peaks == "formulas":  # try to get it from metadata
            try:
                formulas = self.metadata['formula']
            except:
                print("No formula found in metadata. Using numeric labels instead.")
                label_peaks = "numeric"
        elif label_peaks == "numeric": 
            # use rounded mass values as labels
            formulas = [str(round(mass, 2)) for mass in self.peaks.mz]
        else:
            formulas = np.zeros(len(self.peaks.mz)) #TODO get specific formuals by adding option .formulas!
            formulas = [""]*len(self.peaks.mz) #TODO get specific formuals by adding option .formulas!
        if 'run_file' in self.metadata:
            substance = self.metadata['run_file'].split("_")[-1]

        """Plots a first, e.g. the experimental unknown spectra against the a list of other spectra to compare (upside down)"""
        dict_ind = {}
        # check if metadata "compound_bin" is present
        if "compound_bin" in dir(self):
            compound_bin = self.compound_bin
            # get unique compound bins
            unique_compound_bins = np.unique(compound_bin)
            # get indices of unique compound bins
            for ucb in unique_compound_bins:
                dict_ind[str(ucb)] = np.where(compound_bin == ucb)[0]
            # get colors for each unique compound bin
            colors = [plt.cm.tab10(float(i)/np.max([1,float(len(unique_compound_bins)-1)])) for i in range(len(unique_compound_bins))]
            labels_stem = ["comp. #" +  str(ucb) for ucb in unique_compound_bins]
            labels_lines = ["" for ucb in unique_compound_bins]
        else:
            indices = np.arange(len(self.peaks.mz))
            dict_ind['0'] = indices
            colors = [col_stem]
            unique_compound_bins = [0]
            labels_stem = [label_stem[0]]
            labels_lines = [label_stem[1]]

        ax.set_xlabel(r"$m\//z$ (Da)")
        ax.set_ylabel(r"Norm. Intensity (arb. u.)")
        ax.set(ylim=(0, 1.2))

        for i in range(len(unique_compound_bins)):
            ind_sel = dict_ind[str(unique_compound_bins[i])]
            color_plot = colors[i]

            # plot stem
            #ax.plot(self.peaks.mz[ind_sel], self.peaks.intensities[ind_sel])
            markerline, stemline, baseline = plt.stem(self.peaks.mz[ind_sel], self.peaks.intensities[ind_sel], use_line_collection = True, label = labels_stem[i])
            plt.setp(stemline, linewidth=1.0,color=color_plot)
            plt.setp(baseline, linewidth=0.5, color="black")
            plt.setp(markerline, markersize=1.5,color=color_plot)

            # draw assignments
            for d, l, r in zip(self.peaks.mz, self.peaks.intensities, formulas):
                if label_peaks == "formulas":
                    try: 
                        Formula(r)
                    except:
                        continue
                    ax.annotate(r'${}$'.format(formula_to_latex(r)), xy=(d, l), xytext=(2,0),
                                textcoords="offset points", va = "bottom", ha="center")
                elif label_peaks == "numeric":
                    ax.annotate(r'${}$'.format(r), xy=(d, l), xytext=(2,0),
                                textcoords="offset points", va = "bottom", ha="center")
                else:
                    pass
            # draw experimental lines
            if hasattr(self, "meas_mass_in_m_per_z"):
                exp_masses = np.array(self.meas_mass_in_m_per_z)[ind_sel]
            else:
                exp_masses = self.peaks.mz[ind_sel] # it assumes that if measured mass is not in the meatadata it muss be in the peaks (for example because fragments were not assigned)

            for i_em, xc in enumerate(exp_masses): #TODO substitute by mass + uncertainty!
                if i_em==0: 
                    plt.axvline(x=xc, color = color_plot, linewidth=0.2, alpha = 0.5, label=labels_lines[i])
                else:
                    plt.axvline(x=xc, color = color_plot, linewidth=0.2, alpha = 0.5)
        ax.legend()
        #ax.set_title(r'${}$'.format(formula_to_latex(str(formula_from_namestring(Path(path).name)))))
        ax.set_title(r'${}$'.format(str(substance)))
        # draw experimental lines

# needed by StickSpectra class
def calc_disjoint_union_sets(set_list_arb: list) -> list:
    """calculates for a list of sets with potentially nonempty intersection a list of disjoint union sets of unique values

    e.g. [[1, 5], [3, 3, 4], [8,5]] becomes [[1, 5, 8], [3, 4]]

    # arguments:
    :arg set_list_arb: list of arbitrary sets (float or int)
    :return disjoint_union_sets: list of sorted, disjoint union sets of unique values

    """
    a = set_list_arb
    # TODO insert some input testing
    # processed by grouping algo
    if len(a) == 1:
        disjoint_union_sets = [sorted(list(set(a[0])))]
    elif len(a) == 0:
        disjoint_union_sets = []
    else:
        processed = []
        # first elem
        a_first = a[0]
        # second part
        a_rest = a[1:]
        while len(a_rest) >= 0:
            intersection = True
            # as long as an intersection is detected it could be that during for loop intersecting elements are created before the actual processed one
            # when for loop has already passed it -> do it again, unless it is safe that everything is included
            while intersection:
                intersection = False
                for ind, val in enumerate(a_rest):
                    if set(a_first).intersection(set(val)):
                        a_first = set(a_first).union(set(val))
                        a_rest[ind] = []
                        intersection = True
                # remove empty elements of list
                a_rest = [val for val in a_rest if val]
            # the sublist is fully processed and no more overlap -> put it into processed list
            processed.append(sorted(list(a_first)))
            # if not empty list: prepare for the second remaining element
            if a_rest:
                a_first = a_rest[0]
                a_rest = a_rest[1:]
            else:
                break
        disjoint_union_sets = processed
    return disjoint_union_sets

class StickSpectra(ms.Spectrum):
    """Class for stick spectra, which are spectra with discrete peaks and intensities.
        It inherits from the matchms.Spectrum class.
        It is used for theoretical spectra, which are generated by the fragmentation tree or discretized experimental spectra.
    """

    def __init__(
        self,
        mz: np.array,
        intensities: np.array,
        metadata: Optional[dict] = {"precursor_mz": -1},
        metadata_harmonization: bool = True,
        min_mz: Optional[float] = None,
        max_mz: Optional[float] = None,
    ):
        """Initialize StickSpectra class.
        Arguments:
        mz: np.array
            Array of m/z values.
        intensities: np.array
            Array of intensity values.
        metadata: dict
            Dictionary of metadata.
        metadata_harmonization: bool
            If True, harmonize metadata.
        min_mz: float
            Minimum m/z value.
        max_mz: float   
            Maximum m/z value.
        """
        # sort data to avoid super class complaint if not
        _ind_s = np.argsort(mz)

        if all([(mz[i] % 1 == 0) for i in range(len(mz))]):
            val = True
        else:
            val = False
        metadata["is_unity_mz_spec"] = val
        # convert to float (e. g. necessary if unity spectrum and required by matchms)
        # super: inheritance from parent class: initiates calling the init of the match.ms from which it inherites
        super().__init__(
            mz[_ind_s].astype(float),
            intensities[_ind_s].astype(float),
            metadata,
            metadata_harmonization,
        )

        # inherit from AlpinacData the method to plot a stick spectra
        plot_stem = AlpinacData.plot_stem

            # # get fragments name
            # substance = "unknown"
            # if hasattr(self, "formula") and label_peaks == "formulas":
            #         formulas = self.formula
            # elif not hasattr(self, "formula") and label_peaks == "formulas":  # try to get it from metadata
            #     try:
            #         formulas = self.metadata['formula']
            #     except:
            #         print("No formula found in metadata. Using numeric labels instead.")
            #         label_peaks = "numeric"
            # elif label_peaks == "numeric": 
            #     # use rounded mass values as labels
            #     formulas = [str(round(mass, 2)) for mass in self.peaks.mz]
            # else:
            #     formulas = np.zeros(len(self.peaks.mz)) #TODO get specific formuals by adding option .formulas!
            #     formulas = [""]*len(self.peaks.mz) #TODO get specific formuals by adding option .formulas!
            # if 'run_file' in self.metadata:
            #     substance = self.metadata['run_file'].split("_")[-1]

            # """Plots a first, e.g. the experimental unknown spectra against the a list
            # of other spectra to compare (upside down)"""
            # dict_ind = {}
            # # check if metadata "compound_bin" is present
            # if "compound_bin" in dir(self):
            #     compound_bin = self.compound_bin
            #     # get unique compound bins
            #     unique_compound_bins = np.unique(compound_bin)
            #     # get indices of unique compound bins
            #     for ucb in unique_compound_bins:
            #         dict_ind[str(ucb)] = np.where(compound_bin == ucb)[0]
            #     # get colors for each unique compound bin
            #     colors = [plt.cm.tab10(float(i)/np.max([1,float(len(unique_compound_bins)-1)])) for i in range(len(unique_compound_bins))]
            #     labels_stem = ["comp. #" +  str(ucb) for ucb in unique_compound_bins]
            #     labels_lines = ["" for ucb in unique_compound_bins]
            # else:
            #     indices = np.arange(len(self.peaks.mz))
            #     dict_ind['0'] = indices
            #     colors = [col_stem]
            #     unique_compound_bins = [0]
            #     labels_stem = [label_stem[0]]
            #     labels_lines = [label_stem[1]]

            # ax.set_xlabel(r"$m\//z$ (Da)")
            # ax.set_ylabel(r"Norm. Intensity (arb. u.)")
            # ax.set(ylim=(0, 1.2))

            # for i in range(len(unique_compound_bins)):
            #     ind_sel = dict_ind[str(unique_compound_bins[i])]
            #     color_plot = colors[i]

            #     # plot stem
            #     #ax.plot(self.peaks.mz[ind_sel], self.peaks.intensities[ind_sel])
            #     markerline, stemline, baseline = plt.stem(self.peaks.mz[ind_sel], self.peaks.intensities[ind_sel], use_line_collection = True, label = labels_stem[i])
            #     plt.setp(stemline, linewidth=1.0,color=color_plot)
            #     plt.setp(baseline, linewidth=0.5, color="black")
            #     plt.setp(markerline, markersize=1.5,color=color_plot)

            #     # draw assignments
            #     for d, l, r in zip(self.peaks.mz, self.peaks.intensities, formulas):
            #         if label_peaks == "formulas":
            #             try: 
            #                 Formula(r)
            #             except:
            #                 continue
            #             ax.annotate(r'${}$'.format(formula_to_latex(r)), xy=(d, l), xytext=(2,0),
            #                         textcoords="offset points", va = "bottom", ha="center")
            #         elif label_peaks == "numeric":
            #             ax.annotate(r'${}$'.format(r), xy=(d, l), xytext=(2,0),
            #                         textcoords="offset points", va = "bottom", ha="center")
            #         else:
            #             pass
            #     # draw experimental lines
            #     if hasattr(self, "meas_mass_in_m_per_z"):
            #         exp_masses = np.array(self.meas_mass_in_m_per_z)[ind_sel]
            #     else:
            #         exp_masses = self.peaks.mz[ind_sel] # it assumes that if measured mass is not in the meatadata it muss be in the peaks (for example because fragments were not assigned)

            #     for i_em, xc in enumerate(exp_masses): #TODO substitute by mass + uncertainty!
            #         if i_em==0: 
            #             plt.axvline(x=xc, color = color_plot, linewidth=0.2, alpha = 0.5, label=labels_lines[i])
            #         else:
            #             plt.axvline(x=xc, color = color_plot, linewidth=0.2, alpha = 0.5)
            # ax.legend()
            # #ax.set_title(r'${}$'.format(formula_to_latex(str(formula_from_namestring(Path(path).name)))))
            # ax.set_title(r'${}$'.format(str(substance)))
            # # draw experimental lines

# needed by StickSpectra class
def calc_disjoint_union_sets(set_list_arb: list) -> list:
    """calculates for a list of sets with potentially nonempty intersection a list of disjoint union sets of unique values

    e.g. [[1, 5], [3, 3, 4], [8,5]] becomes [[1, 5, 8], [3, 4]]

    # arguments:
    :arg set_list_arb: list of arbitrary sets (float or int)
    :return disjoint_union_sets: list of sorted, disjoint union sets of unique values

    """
    a = set_list_arb
    # TODO insert some input testing
    # processed by grouping algo
    if len(a) == 1:
        disjoint_union_sets = [sorted(list(set(a[0])))]
    elif len(a) == 0:
        disjoint_union_sets = []
    else:
        processed = []
        # first elem
        a_first = a[0]
        # second part
        a_rest = a[1:]
        while len(a_rest) >= 0:
            intersection = True
            # as long as an intersection is detected it could be that during for loop intersecting elements are created before the actual processed one
            # when for loop has already passed it -> do it again, unless it is safe that everything is included
            while intersection:
                intersection = False
                for ind, val in enumerate(a_rest):
                    if set(a_first).intersection(set(val)):
                        a_first = set(a_first).union(set(val))
                        a_rest[ind] = []
                        intersection = True
                # remove empty elements of list
                a_rest = [val for val in a_rest if val]
            # the sublist is fully processed and no more overlap -> put it into processed list
            processed.append(sorted(list(a_first)))
            # if not empty list: prepare for the second remaining element
            if a_rest:
                a_first = a_rest[0]
                a_rest = a_rest[1:]
            else:
                break
        disjoint_union_sets = processed
    return disjoint_union_sets

class StickSpectra(ms.Spectrum):
    """Class for stick spectra, which are spectra with discrete peaks and intensities.
        It inherits from the matchms.Spectrum class.
        It is used for theoretical spectra, which are generated by the fragmentation tree or discretized experimental spectra.
    """

    def __init__(
        self,
        mz: np.array,
        intensities: np.array,
        metadata: Optional[dict] = {"precursor_mz": -1},
        metadata_harmonization: bool = True,
        min_mz: Optional[float] = None,
        max_mz: Optional[float] = None,
    ):
        """Initialize StickSpectra class.
        Arguments:
        mz: np.array
            Array of m/z values.
        intensities: np.array
            Array of intensity values.
        metadata: dict
            Dictionary of metadata.
        metadata_harmonization: bool
            If True, harmonize metadata.
        min_mz: float
            Minimum m/z value.
        max_mz: float   
            Maximum m/z value.
        """
        # sort data to avoid super class complaint if not
        _ind_s = np.argsort(mz)

        if all([(mz[i] % 1 == 0) for i in range(len(mz))]):
            val = True
        else:
            val = False
        metadata["is_unity_mz_spec"] = val
        # convert to float (e. g. necessary if unity spectrum and required by matchms)
        # super: inheritance from parent class: initiates calling the init of the match.ms from which it inherites
        super().__init__(
            mz[_ind_s].astype(float),
            intensities[_ind_s].astype(float),
            metadata,
            metadata_harmonization,
        )

    # inherit from AlpinacData the method to plot a stick spectra
    plot_stem = AlpinacData.plot_stem

    @classmethod
    def _get_property(cls, _property_name):
        """get properties from metadata attribute if not already setted."""
        if (
            getattr(cls, _property_name) is None
        ):  # if internal property uncertainty not set yet, it will set it
            if _property_name in cls:
                property = cls._property_name
                setattr(cls, name=_property_name, value=property)
            raise ValueError(f"{_property_name} not set in {cls}.")
        return getattr(cls, _property_name)

    # @_set_property_to_metadata.setter # sets properties, so that x.property = y works
    @classmethod
    def _set_property(cls, _property_name, _value):
        cls._property_name = _value
        setattr(cls, _property_name, _value)

    def bin_around_center_position_of_ref_spectra(
        self, theo_spec: StickSpectra, tollevel: float
    ) -> tuple[np.array, np.array]:
        """function generates bins around reference peaks (likely theoretical) and assigns peaks of measured spectrum to those mass bins
        if they are within tolerance level

        # arguments:
        :arg theo_spec: spectrum used as reference: on all peaks of the spectra a bin of a width 2*tollevel is created
        :arg tollevel: half width of tolerance delta mz

        # returns:
        :return bin_theo_eval: array of bin, center mass and intensity of reference data (if several bins overlap these bins are merged into a broader one)
        :return bin_meas_eval: array of bin, center mass and intensity of input data

        """
        theo_spec_m = theo_spec.peaks.mz
        theo_spec_I = theo_spec.peaks.intensities
        meas_spec_m = getattr(getattr(self, "peaks"), "mz")
        meas_spec_I = getattr(getattr(self, "peaks"), "intensities")
        # the bins are created using the theoretical spectrum/ spectrum1 (where the peaks are).
        bin_l = theo_spec_m - tollevel
        bin_r = theo_spec_m + tollevel
        bin_h = [[bin_l[i], bin_r[i]] for i in range(len(bin_l))]
        # if two bins overlap, combine these bins to one bigger one
        for i in range(len(bin_l) - 1):
            # if right border of lighter mass interval is larger than left border of next heavy mass interval: combin
            if bin_h[i][1] > bin_h[i + 1][0]:
                # enlarge left border of next interval
                bin_h[i + 1][0] = bin_h[i][0]
                # remove current interval
                bin_h[i][0] = None
                bin_h[i][1] = None
        bin_hh = [i for i in bin_h if i[0]]
        # flatten intervall
        bin = list(chain(*bin_hh))
        # getting indici of values in the corresponding intervals: 0 is before, odd are between of mass peaks (no signal), even are in between
        ind_bin_theo = np.digitize(theo_spec_m, bin)
        ind_bin_meas = np.digitize(meas_spec_m, bin)
        # get all elements witch are within the theoretical peaks, other peaks in the experimental spectrum are ignored
        # in an experimental spectrum with several independent (no matching peaks within an interval of same width) spectra
        # scaling can be performed independently. Be careful, if mass peaks of two independent spectra are lying in the same interval, adapt tollevel.
        bin_theo_eval = [
            list([ind_bin_theo[i], theo_spec_m[i], theo_spec_I[i]])
            for i in range(len(theo_spec_m))
            if ind_bin_theo[i] % 2 != 0
        ]
        bin_meas_eval = [
            list([ind_bin_meas[i], meas_spec_m[i], meas_spec_I[i]])
            for i in range(len(meas_spec_m))
            if ind_bin_meas[i] % 2 != 0 and ind_bin_meas[i] in ind_bin_theo
        ]
        return np.array(bin_theo_eval), np.array(bin_meas_eval)

    def scale_fac_rel_to(self, ref_spectrum: StickSpectra, tollevel: float) -> float:
        """determines one global scaling factor for the input stickspectra relative to a reference spectrum by least-square optimization
        of the intensity difference of input and ref spectra within the bins of the reference spectrum

        # arguments:
        :arg ref_spec: spectrum used as reference: only intensities at bins located at the ref_spec peak positions are taken into account.
        :arg tollevel: half width of tolerance delta mz

        # returns:
        :return scale_opt: optimized global scaling parameter

        """
        # generate intervals/bins
        bin_theo_eval, bin_meas_eval = self.bin_around_center_position_of_ref_spectra(
            theo_spec=ref_spectrum, tollevel=tollevel
        )
        # get scaling parameter by optimization:
        no_unique_bins = len(
            np.unique([bin_theo_eval[i][0] for i in range(len(bin_theo_eval))])
        )
        # vectorize intensities
        int_vec_theo = np.zeros(no_unique_bins)
        int_vec_meas = np.zeros(no_unique_bins)
        # sum if several in one bin and form a vector, each position corresponding to a theoretical bin.
        for i in range(len(bin_theo_eval)):  # i=0
            int_vec_theo[int((bin_theo_eval[i][0] - 1) / 2)] += bin_theo_eval[i][2]
        for i in range(len(bin_meas_eval)):  # i=0
            int_vec_meas[int((bin_meas_eval[i][0] - 1) / 2)] += bin_meas_eval[i][2]
        # max scaling of signals in same bin
        # remove zeros
        mask = int_vec_theo != 0
        int_vec_theo = int_vec_theo[mask]
        int_vec_meas = int_vec_meas[mask]
        scale_est = max(np.array(int_vec_meas) / np.array(int_vec_theo))
        param = Parameters()
        param.add("scale", value=scale_est)
        # optimization: minimizing cost function to get relative values
        # last variable is 1 - all others.
        minner = Minimizer(
            StickSpectra._cost_scale, param, fcn_args=(int_vec_theo, int_vec_meas)
        )
        # optimize by using Limited-memory BFGS,
        # an algorithm determining the optimum based on the derivative of N-dimensional potential surface
        peak_optimized = minner.minimize(method="lbfbsg")
        scale_opt = peak_optimized.params["scale"].value
        return scale_opt

    def vectorize_intensities(self, ref_spec, tollevel):
        """vectorizes intensities in bins with width = 2*tollevel defined by ref_spec and returns intenisity vectors for reference and
            input spectra (self)

        # arguments:
        :arg ref_spec: spectrum used as reference: only intensities at bins located at the ref_spec peak positions are taken into account
        :arg tollevel: half width of tolerance delta mz

        # returns:
        :return int_vec_ref: intensity vector of reference spectrum
        :return int_vec_theo: intensity vector of input spectrum (self)
        """

        # generate intervals/bins
        bin_ref_eval, bin_self_eval = self.bin_around_center_position_of_ref_spectra(
            theo_spec=ref_spec, tollevel=tollevel
        )
        # get scaling parameter by optimization:
        no_unique_bins = len(
            np.unique([bin_ref_eval[i][0] for i in range(len(bin_ref_eval))])
        )
        # vectorize intensities
        int_vec_ref = np.zeros(no_unique_bins)
        int_vec_self = np.zeros(no_unique_bins)
        # sum if several in one bin and form a vector, each position corresponding to a theoretical bin.
        for i in range(len(bin_ref_eval)):  # i=0
            int_vec_ref[int((bin_ref_eval[i][0] - 1) / 2)] += bin_ref_eval[i][2]
        for i in range(len(bin_self_eval)):  # i=0
            int_vec_self[int((bin_self_eval[i][0] - 1) / 2)] += bin_self_eval[i][2]
        return int_vec_ref, int_vec_self

    def fit_absolute_stick_spectra_contributions(
        self, ref_spectra: list[StickSpectra] | StickSpectra, tollevel: float
    ) -> float:
        """determines one global scaling factor for the input stickspectra relative to a reference spectrum by least-square optimization
        of the intensity difference of input and ref spectra within the bins of the reference spectrum

        # arguments:
        :arg ref_spec: spectrum used as reference: only intensities at bins located at the ref_spec peak positions are taken into account.
        :arg tollevel: half width of tolerance delta mz

        # returns:
        :return scale_opt: optimized global scaling parameter

        """
        if not isinstance(ref_spectra, list):
            scale = [StickSpectra.scale_fac_rel_to(self, ref_spectra, tollevel)]
            # print("scale")
            # print(scale)
            return scale
        else:
            # normalize (for fitting all parameters should have same order of magnitude, if weighted by mz distance (normalized) a normalization is also helpful)
            # will lead to fitting parameters between 0 and 1
            norm_fac = max(self.peaks.intensities)
            meas_spec = StickSpectra(self.peaks.mz, self.peaks.intensities / norm_fac)
            # for more than one theoretical spectra: do superposition
            # create as reference spectrum a superposition of all measured and theoretical spectra
            # that will define all mz bins considered for building the intensity vectors (a vector of intensities where each entry correspond to
            # intensity I(i) in bin mz_interval(i)
            ref_superpos = meas_spec.superpose(spectra=ref_spectra)
            # vectorize all meas. and theor. spectra relative to that ref. spectra:
            int_vec_ref, int_vec_meas = meas_spec.vectorize_intensities(
                ref_superpos, tollevel
            )
            # get single contributions for each ref_spectrum of ref_spectra (e.g. isotopic spectra)
            all_ref_intensities_vec = []
            all_mean_delta_m = np.empty(shape=(0, len(ref_spectra)))

            for i, spec in enumerate(ref_spectra):
                # print(spec)
                # vectorize all intensities resp. to common mz basis of meas. and ref. spectra
                int_vec_ref_spec = spec.vectorize_intensities(ref_superpos, tollevel)[1]
                all_ref_intensities_vec.append(int_vec_ref_spec)
                # calculate average deviation of theo. and meas. peaks
                if not meas_spec.get("is_unity_mz_spec"):
                    (
                        diff,
                        weight_diff,
                        spectrum_meas_calib,
                        delta_mass_opt,
                    ) = meas_spec.get_delta_m_for_given_reference_spectrum(
                        spec_theo_ref=spec, tollevel=tollevel
                    )
                    all_mean_delta_m = np.append(
                        all_mean_delta_m, np.average(diff, weights=weight_diff)
                    )
                else:
                    all_mean_delta_m = np.append(all_mean_delta_m, 1)

                # normalize
                all_mean_delta_m_norm = abs(
                    all_mean_delta_m * 1 / np.max(abs(all_mean_delta_m))
                )

            # get initial input parameter as maximal contributions to measured spectramax scaling of signals in same bin
            max_contrib = np.empty(shape=(len(all_ref_intensities_vec)))
            for ind, int_vec in enumerate(all_ref_intensities_vec):
                # remove zeros
                mask_nonzero1 = int_vec != 0 # mask of non zero elements
                mask_nonzero2 = int_vec_meas != 0 # mask of non zero elements
                int_vec_non_zero = int_vec[mask_nonzero1*mask_nonzero2]
                int_vec_meas_non_zero = int_vec_meas[mask_nonzero1*mask_nonzero2]
                if len(int_vec_non_zero) == 0:
                    max_contrib[ind] = 0
                else:
                    ratio = int_vec_meas_non_zero / int_vec_non_zero
                    max_contrib[ind] = ratio[np.isfinite(ratio)].min()
                # print(ratio)
            # in case only a few spectra contribute which do not cover the full mass range (and are maybe even fitted in a seperate step)
            # set measured contribution which can not be covered by isotope spectra equal to zero:
            # sum to detect 0 equiv to only measured spectra
            ref_only_superpos = np.sum(all_ref_intensities_vec, axis =0)
            ind_where_theo_peak = np.nonzero(ref_only_superpos>0)
            int_vec_meas_with_theo_overlap = int_vec_meas[ind_where_theo_peak]
            all_ref_intensities_vec_same_as_meas = [ref[ind_where_theo_peak] for ref in all_ref_intensities_vec]

            params = Parameters()
            for i in range(len(ref_spectra)):
                # parameters per peak:
                # max values is 5% larger plus might have a peak below noise
                #params.add('k'+str(i), value=1/norm_fac*float(ref_spectra[i].get("iso_idx_list_k_initial")), min = 0, max = 1)
                params.add("k" + str(i), value=max_contrib[i], min=0, max=1)
            # print(params)
            # optimization: minimizing cost function to get relative values
            # last variable is 1 - all others.
            minner = Minimizer(
                StickSpectra._cost_f_iso_discrete_multi,
                params,
                #fcn_args=(all_ref_intensities_vec, int_vec_meas, all_mean_delta_m_norm),
                fcn_args=(all_ref_intensities_vec_same_as_meas, int_vec_meas_with_theo_overlap, all_mean_delta_m_norm),
            )
            # optimize by using Limited-memory BFGS,
            # an algorithm determining the optimum based on the derivative of N-dimensional potential surface
            peak_optimized = minner.minimize(method="lsbgfs")
            # generate list of optimized parameters
            k_opt = np.empty(shape=(0, len(ref_spectra)))
            for i in range(len(ref_spectra)):
                k_opt = np.append(k_opt, peak_optimized.params["k" + str(i)].value)
            k_opt_abs = [k * norm_fac for k in k_opt]
            return k_opt_abs

    def in_mz_range(self, mz_min, mz_max):  # cls = spectrum_meas_full
        """determines stick spectrum within m/z borders and returns it

        # arguments:
        :mz_min: m/z lower limit
        :mz_max: m/z upper limit
        # returns:
        :return spec_in_range: stick spectrum within borders
        :return ind_perm: index of peaks within borders
        """
        mz_peak = getattr(getattr(self, "peaks"), "mz")
        mz_int = getattr(getattr(self, "peaks"), "intensities")
        # index of elements within border
        ind_perm_h = np.argwhere((mz_min <= mz_peak) * ((mz_max >= mz_peak)))
        ind_perm = list(chain(*ind_perm_h))
        spec_in_range = StickSpectra(mz=mz_peak[ind_perm], intensities=mz_int[ind_perm])
        return spec_in_range, ind_perm

    @classmethod
    def linear_dependence_of_spectra(
        cls,
        list_spectra: list,
        tolerance_level: float = 0.005,
        threshold_score_for_lin_dep: float = 0,
    ):
        """calculates for a list of spectra the inici of linear dependent and independent spectra

        Decides on the base of a minimal score threshold if spectra depend linearly on each others
        using cosine similarity (matchms functionality)

        # arguments:
        :arg list_spectra: list of spectra in matchms class formate
        :arg threshold_score_for_lin_dep: minimum matchms score to decide whether spectra are linear dependent
        :arg tolerance level: m/z tolerance for which m/z values are attributed to the same peak
        # returns:
        :return lin_dep_spectra: list of list of linear dependent spectra (index);
        :return lin_indep_spectra: list of linear independent spectra (index);
        """
        if len(list_spectra) == 1:
            lin_dep_spectra = [[], []]
            lin_indep_spectra = [[list_spectra[0]], [0]]
            return lin_dep_spectra, lin_indep_spectra
        else:
            # tolerance_level = tollevel
            similarity_measure = CosineGreedy(tolerance=tolerance_level)
            # scores 2D:
            scores2D = calculate_scores(
                references=list_spectra,
                queries=list_spectra,
                similarity_function=similarity_measure,
                is_symmetric=True,
            )
            # score matrix to get the correlation
            score_matr = scores2D.scores.to_array('CosineGreedy_score')
            logger.debug(f"{score_matr=}")
            # get triangle above the diagonal
            offdiag_upper = np.triu(score_matr, k=1)
            # get row,col indici of all elements with values above a threshold or non zero: zero means no overlap
            all_spectra = list(range(score_matr.shape[0]))
            # above threshold
            all_entries_above_tresh = np.argwhere(
                offdiag_upper > threshold_score_for_lin_dep
            )
            # convert them to sets
            spec_lin_dep = [set(x) for x in all_entries_above_tresh.tolist()]
            # remove empy list from list
            # spectra which must be fitted as a superposition
            ind_dep_spectra = calc_disjoint_union_sets(spec_lin_dep)
            # spectra which can be fitted independently (difference set of dependent and
            ind_indep_spectra = list(
                set(all_spectra).difference(
                    set(list(chain(*ind_dep_spectra)))
                )
            )
            # ind is grouped, e.g. [[1,2,3],[4,8]] -> also group sepctrum
            spec_gr = []
            for ii in range(len(ind_dep_spectra)):
                spec_gr.append([list_spectra[i] for i in ind_dep_spectra[ii]])
            lin_dep_spectra_gr = [spec_gr, ind_dep_spectra]
            lin_indep_spectra = [
                [list_spectra[i] for i in ind_indep_spectra],
                ind_indep_spectra,
            ]
            return lin_dep_spectra_gr, lin_indep_spectra

    @classmethod
    def aggregate_np(cls, mz: np.array, int: np.array, tolerance: float = 0, return_mean_for_multiple_in_tolerance: bool = True):
        ind_sort = np.argsort(mz)
        mz_sorted = mz[ind_sort]
        int_sorted = int[ind_sort]
        cummulated_int = np.empty(shape=[0, 1])
        unique_mz = np.unique(mz_sorted)
        if return_mean_for_multiple_in_tolerance: # option since 13.06.23: if several peaks are too close (within tolerance) take mean of them 
            mz_list_sets = []
            for uniq_val in unique_mz:
                # get intersection of all mz values within tolerance
                ind_intersect = np.logical_and(unique_mz>= uniq_val - tolerance, unique_mz <= uniq_val + tolerance)
                # get mean of all mz values with ind_intersect == True
                mz_list_sets.append(np.mean(unique_mz[ind_intersect]))  
            unique_mz = np.unique(np.array(mz_list_sets))

        for uniq_val in unique_mz:
            all_int = [
                int_sorted[ind]
                for ind, val in enumerate(mz_sorted)
                if (
                    mz_sorted[ind] >= uniq_val - tolerance
                    and mz_sorted[ind] <= uniq_val + tolerance
                )
            ]
            cummulated_int = np.append(cummulated_int, sum(all_int))
        return unique_mz, cummulated_int

    def aggregate_by_mass_or_mass_interval(self, tolerance: float = 0):
        """aggregates a spectrum of class StickSpectra if multiple peaks of same m/z or in same bin (mz +- tollerance) are present

        # arguments:
        :arg tolerance: m/z tolerance for which m/z values are attributed to the same m/z
        # returns:
        :return spec_sorted_aggregated: aggregated stick spectrum
        """
        # acummulate/sum elements of same mass
        mz_sorted = self.peaks.mz
        int_sorted = self.peaks.intensities
        unique_mz, cummulated_int = StickSpectra.aggregate_np(
            StickSpectra, mz_sorted, int_sorted, tolerance
        )
        spec_sorted_aggregated = StickSpectra(mz=unique_mz, intensities=cummulated_int)
        return spec_sorted_aggregated

    #
    def superpose(
        self,
        spectra: StickSpectra | list[StickSpectra],
        weight_self: float | int = 1,
        weight_spectra: list[float] | float | list[int] | int = 1,
        tollevel_aggregate: float = 0,
    ) -> StickSpectra:
        """superpose a spectrum of class StickSpectra with one or more other stick spectra

        # arguments:
        :arg spectra: (list of) spectra/um which should be superimposed
        :arg weight_self:  weight with which the contributing intensities of input/self
        :arg weight_spectra:  (list of) weight/s with which the contributing intensities of argument spectra
        :arg tollevel_aggregate: m/z tolerance for which m/z values are attributed to the same m/z while aggrediation
        # returns:
        :return superposition: superposition stick spectrum
        """
        # set to same dimensions if only single spectrum is given
        if not isinstance(spectra, list):
            spectra = [spectra]
        # number of spectra for superposition
        n_spectra = len(spectra)
        # if only one value for weight is given, fill up array with this value, default one
        if isinstance(weight_spectra, float | int):
            weightn = np.empty(n_spectra)
            weightn.fill(float(weight_spectra))
            weight_spectra = weightn
        elif len(weight_spectra) != n_spectra:
            raise ValueError(
                f"weights array length {len(weight_spectra)} differ from spectra array length {n_spectra}."
            )
        # if weight = 0 it should not add the zero
        if weight_self != 0:
            all_mz = self.peaks.mz
            all_int = self.peaks.intensities * weight_self
        else:
            all_mz = np.empty(shape=[0, 1])
            all_int = np.empty(shape=[0, 1])

        # accumulate peaks
        for i in range(n_spectra):
            if weight_spectra[i] != 0:
                all_mz = np.append(all_mz, spectra[i].peaks.mz)
                all_int = np.append(
                    all_int, spectra[i].peaks.intensities * weight_spectra[i]
                )
        # it is assumed you pass at least one valid spectrum with a weight > 0
        superpos_raw_mz, superpos_raw_int = StickSpectra.aggregate_np(
            mz=all_mz, int=all_int, tolerance=tollevel_aggregate
        )
        # cummulate intensities of same mass (interval)
        # superposition = superpos_raw.aggregate_by_mass_or_mass_interval(
        #    tolerance=tollevel_aggregate
        # )
        superposition = StickSpectra(mz=superpos_raw_mz, intensities=superpos_raw_int)
        return superposition

    def matching_score_with(
        self, spectrum_to_compare: ms.Spectrum, tolerance_level: float = 0.005
    ) -> float:
        """calculates the matching (number of commun peaks and cosine greedy
        scores within a tolerance mass shift (shift in m/z)) and returns total score
        # arguments:
        :arg second_spectrum: spectrum in matchms class formate for comparison
        :arg tolerance level: m/z tolerance for which m/z values are attributed to the same peak
        # returns:
        :return score: cosine similarity score
        """
        # strategy: reduce dimension of fitting problem by checking which spectra depend from each other
        # to do so do a 2D match ms search and group into groups dependent of each others
        # the sub group will form superpositions of Sum(k_i*spectrum_i) and an do an optimization on them.
        # fit seperately for a) superposition of linear dependent spectra and b) independent spectra

        similarity_measure = CosineGreedy(tolerance=tolerance_level)
        scores0D = calculate_scores(
            references=[self],
            queries=[spectrum_to_compare],
            similarity_function=similarity_measure,
            is_symmetric=False,
        )
        score = scores0D.scores[0][0]["score"]
        # return -(score-1) # if 0 matching -> 1; if 1 mathing -> 0
        return score

    def fit_relative_stick_spectra_contributions(
        self,
        iso_stick_spectra: list[StickSpectra],
        tollevel: float,
        iso_idx_list_k_initial: list,
    ):
        """fits/optimizes the relative contributions (normalized to 1)
        of all spectra given by iso_stick_spectra to a spectra of class StickSpectra.
        # arguments:
        :arg iso_stick_spectra: list of spectra contributing (potentially) to input spectrum
        :arg tolerance level: m/z tolerance for which m/z values are attributed to the same peak
        # returns:
        :return k_opt: optimized relative contributions
        :return superpos_spectrum: resulting optimized superposition spectrum
        """
        # TODO remove this function, not needed anymore
        no_iso_spectra = len(iso_stick_spectra)
        # initialize fit parameters
        params = Parameters()
        # estimates for all peaks the minimum value of k and adds it to parameter list, as only relative behaviour is optimized last one is 1 - sum(k_i(rest))
        for idx_peak in range(len(iso_stick_spectra) - 1):
            # parameters per peak:
            params.add(
                "k" + str(idx_peak),
                value=iso_idx_list_k_initial[idx_peak] / sum(iso_idx_list_k_initial),
                min=0,
                max=1,
            )

        # optimization: minimizing cost function to get relative values
        # last variable is 1 - all others.
        # Todo kaho: remove hardcoded 3! HARD-CODED
        minner = Minimizer(
            StickSpectra._cost_f_iso_discrete,
            params,
            fcn_args=(iso_stick_spectra, self, tollevel),
        )
        # the cost value is the matchms score -1 it will not lead to reliable uncertainties in least square evaluation
        # implemented as it is, it is a fit with one data point ...
        # peak_optimised = minner.minimize(method = 'leastsq') -> not good, if to little data points available
        peak_optimised = minner.minimize(method="lbfbsg")

        # # optimized intensities
        # superpose on last isotopic profile the other profiles
        k_opt_all_but_last_spectra = np.array(
            [
                peak_optimised.params["k" + str(i)].value
                for i in range(no_iso_spectra)
                if i < no_iso_spectra - 1
            ]
        )
        k_opt_last_spectra = 1 - sum(k_opt_all_but_last_spectra)
        k_opt = np.append(k_opt_all_but_last_spectra, k_opt_last_spectra)
        # iso_stick_spectra_last_weighted = StickSpectra(
        #     mz=iso_stick_spectra[-1].peaks.mz,
        #     intensities=iso_stick_spectra[-1].peaks.intensities * k_opt_last_spectra,
        # )
        # TODO kaho HARD-CODED
        superpos_spectrum = iso_stick_spectra[-1].superpose(
            iso_stick_spectra[:-1],
            weight_self=k_opt_last_spectra,
            weight_spectra=k_opt_all_but_last_spectra,
            tollevel_aggregate=0,
        )
        return k_opt, superpos_spectrum

    def get_delta_m_for_given_reference_spectrum(self, spec_theo_ref, tollevel):
        (
            bin_theo_eval,
            bin_meas_eval,
        ) = self.bin_around_center_position_of_ref_spectra(
            theo_spec=spec_theo_ref,
            tollevel=tollevel,
        )
        # if at least one overlapping bin of theo and meas. data is found
        diff = []
        weight_diff = []
        spectrum_meas_calib = []
        delta_mass_opt = np.NaN

        if len(bin_meas_eval) > 0:
            for ii in range(len(bin_meas_eval)):
                for jj in range(len(bin_theo_eval)):
                    if bin_meas_eval[ii][0] == bin_theo_eval[jj][0]:
                        diff.append(bin_meas_eval[ii][1] - bin_theo_eval[jj][1])
                        # if the measured or theoretical intensity is heigher, the difference weights stronger
                        # low measured mass peaks were liekly not fitted well in pre-processing
                        # low theoretical peaks could also likely belong to a measured peak not appearing in the measured spectra and are attributed to a wrong peak
                        weight_diff.append(bin_meas_eval[ii][2] * bin_theo_eval[jj][2])
                        # should be also defined in the rest of the code as rel. to theo/true, if positive: mass_meas too high
            delta_mass_opt = np.average(diff, weights=weight_diff)
            # Note: if not well calibrated data, in a while loop, the measured spectra could be corrected and the tolerance interval could be reduced and again optimization
            spectrum_meas_calib = StickSpectra(
                mz=self.peaks.mz - float(delta_mass_opt),
                intensities=self.peaks.mz,
            )
        else:
            pass
            # TODO set logger
            print(
                f"No overlap of theoretical {bin_theo_eval} and measured peaks {bin_meas_eval} is detected. Acceptance interval is extended by factor 1.5"
            )
        return diff, weight_diff, spectrum_meas_calib, delta_mass_opt

    @classmethod
    def _cost_f_iso_discrete(
        cls,
        params: Parameters,
        iso_stick_spectra: StickSpectra | list[StickSpectra],
        spectrum_meas: StickSpectra,
        tollevel: float,
    ):
        """cost function returning the matching of the input spectrum and
        the superposition of spectra (iso_stick_spectra) weighted by an optimization input parameter"""
        # number of isotopologue mass spectra belonging to the node (could go outside cost)
        number_of_stick_spectra_eval = len(iso_stick_spectra)
        # parameter list
        paramlist = [
            params["k" + str(i)] for i in range(number_of_stick_spectra_eval - 1)
        ]

        last_spectrum = iso_stick_spectra[number_of_stick_spectra_eval - 1]
        spectra_to_superpose = iso_stick_spectra[: number_of_stick_spectra_eval - 1]
        # superposition without aggregating
        # last weight/ k is what is left from 100%, weight is 1 - summed percentage of all other contributions
        # other weights/ ks are the fitting parameters directly
        spectrum_model_ms = last_spectrum.superpose(
            spectra=spectra_to_superpose,
            weight_self=(1 - sum(paramlist)),
            weight_spectra=np.array(paramlist),
            tollevel_aggregate=tollevel,
        )
        # also consider delta m
        if not spectrum_meas.get("is_unity_mz_spec"):
            (
                diff,
                weight_diff,
                spectrum_meas_calib,
                delta_mass_opt,
            ) = spectrum_meas.get_delta_m_for_given_reference_spectrum(
                spec_theo_ref=spectrum_model_ms, tollevel=tollevel
            )
            # high value if masses far away from the center have a high intensity, masses next to center will be preferred
            # divided by tollevel to make sure that is around 1
            delta_m_weighted_sum = (
                sum(abs(np.array(diff) * np.array(weight_diff))) / tollevel
            )
        else:
            delta_m_weighted_sum = 1
        cost = np.sqrt(
            abs(
                spectrum_meas.matching_score_with(
                    spectrum_model_ms, tolerance_level=tollevel
                )
                - 1
            )
        )
        return cost + (delta_m_weighted_sum)

    @classmethod
    def _cost_scale(cls, param, int_vec_theo, int_vec_meas):
        """cost function returning the absolute difference value between two vectors,
        where each position refers to the same mass bin: the static (, likely measured) spectrum and a theoretical spectrum
        and therefore the absolute contribution of the theoretical spectrum"""
        diff = abs(int_vec_meas - param["scale"] * int_vec_theo)
        return diff

    @classmethod
    def _cost_f_iso_discrete_multi(
        cls,
        params,
        all_ref_intensities_vec: np.array,
        int_vec_meas: np.array,
        all_mean_delta_m_norm: np.array,
    ):
        """
        cost is deviation from all intensities in same mass bin multiplied by deviation of
        average mass deviation from 0 (normalized to 1) multiplied by a weight
        parameters are the relative contributions of the isotopic spectra
        # it optimizes the relative contributions of multiple spectra contributions (a global scaling factor must be found with _cost_scale)
        """

        int_all_ref = params["k0"] * all_ref_intensities_vec[0]
        all_weight = all_mean_delta_m_norm[0] * params["k0"]
        #        if params["k0"] > 1E-4:
        #            w2 = 1
        #        else:
        #            w2 = 0
        for i in range(len(all_ref_intensities_vec))[1:]:
            # print(int_all_ref)
            int_all_ref = (
                int_all_ref + params["k" + str(i)] * all_ref_intensities_vec[i]
            )
            all_weight = all_weight + all_mean_delta_m_norm[i] * params["k" + str(i)]
        #            if params['k'+str(i)] > 1E-4:
        #                w2+=1
        # fitting failes if one vec is zero
        norm = np.linalg.norm(int_vec_meas) * np.linalg.norm(int_all_ref)
        if norm == 0:
            cos_sim = 1
        else:
            cos_sim = 1 - np.dot(int_vec_meas, int_all_ref) / (norm)
        diff = np.abs(int_vec_meas - int_all_ref)
        # weight difference in intensity matching 80% and total deviation of mass 20%
        cost = (
            0.45 * cos_sim + 0.45 * diff + 0.035 * all_weight
        )  # +0.015*w2/len(all_ref_intensities_vec)
        return cost

def formula_from_namestring(namestring: str):
    """ return formula from a namestring, as used in ALPINAC input/output files """
    name_parts = re.split(r'[., \-!?:_]+', namestring)
    name = None
    for piece in name_parts:
        try: 
            name = Formula(piece)
        except:
            continue
    return name

def alpinac_input_2_mgf_converter(alpdata:AlpinacData.input)->ms.spectrum: #TODO: into AlpinacData class
    """ under development: should be a method which converts AlpinacData.input Object into mfg"""
    #TODO: next convert all files from 
    txt_path = r"C:\Users\kaho\polybox\ALPINAC\test_kaho\test_data\test_minimal_mfg.txt"
    spectrums = list(load_from_mgf(txt_path))
    type(spectrums)
    spectrums[0].metadata

# for visualization of pyside
def get_list_of_het_atoms_from_formula(formula:str):
    """ returns a list of heterogenous atoms from a chemical formula as default input parameter for ALPINAC."""
    el_list = [str(el) for el in elements] # to init file? requests module periodictabl;e
    hetatoms = [el for el in el_list if el in formula]
    #check if e.g. H was confused with He by comparing if 1st character from two elements occurs at same pos
    hetatoms.sort(key=lambda s: len(s), reverse=True)
    occurence = []
    for el in hetatoms:# el = hetatoms[0]
        #position on which element is found
        occurence_el = [m.start() for m in re.finditer(el, formula)] #el = 'C' 
        # if at the same position already a 2-character elements is located, the following 1 char-element is removed Ca -> 0 later C-> 0, remove C
        if set(occurence_el).issubset(set(occurence)):  
            hetatoms.remove(el)
        else:
            occurence = occurence + occurence_el
    return hetatoms

def smiles_to_formula(smiles:str):
    """ converts smiles code into chemical formula"""
    mol = Chem.MolFromSmiles(smiles)
    formula = CalcMolFormula(mol)
    return formula

def formula_to_latex(formula:str):
    """ converts a chemical formula into latex code """
    latex_formula = Formula(formula).latex
    return latex_formula

def get_list_of_het_atoms_from_smiles(smiles:str):
    """ returns list as default input parameter for ALPINAC"""
    mol = Chem.MolFromSmiles(smiles)
    hetatoms = {atom.GetSymbol() for atom in mol.GetAtoms()} #TODO
    return hetatoms

def draw_chemical_struct(ax, spectrum: ms.Spectrum):
    """draws an image of a spectrum with a smiles code into ax"""
    smiles = spectrum.get("smiles")
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        formula = smiles_to_formula(smiles)
        latex_formula = formula_to_latex(formula)
        # CalcMolFormula(m)
        pil_img = Draw.MolToImage(mol)
        plt.imshow(pil_img)
        plt.text(x=0.1, y=0.9, s=f'${latex_formula}$',transform=plt.gca().transAxes)
    else: return None
    
def halocarbname2formula(name:str, except_dict: dict = {})->str:
    """ converts a halocarbon name into a chemical formula
    # arguments:
    :arg name: halocarbon name
    :arg except_dict: dictionary with exceptions
    # returns:
    :return formula: chemical formula
    """
    # TODO figure out HCFO name convention
    name_PubChem = [] #a structure information can be obtained by the letters behind
    #There are some substances which do not fit to conventions, e.g. if more than 9 atoms -> 10 but functions does not know if HFC5111 is HFC5-11-1 or HFC5-1-11 and Myriam was not using -: so use exception dict
    [no_carbs, no_hydrs, no_chlorines, no_fluors, no_bromines] = [0,0,0,0,0]
    # homogenise name
    name = name.replace("-",'').replace("-","").replace(" ","")
    # replace small letters at the end of expression if the follow a hyphen followed by any number, e.g. HFC-1234yf -> HFC-1234
    name = re.sub(r'(?<=\d)[a-z]+$', '', name)
    if name in except_dict.keys():
        [formula, name_PubChem] = [except_dict[name], ""]
    # extrac chemical group
    elif name.lower().startswith("halon"):
        # halons (fluorocarbons with Br)
        group = "halon"
        name = name.lower()
        atom_info = name.split(group)[1]
        def_num = "".join([char for char in atom_info if char.isdigit()])
        [no_carbs, no_fluors, no_chlorines, no_bromines] = [int(i) for i in def_num]
        name_PubChem = name.replace("halon","Halon-")
    elif name.lower().startswith("hfo"):
        # halons (fluorocarbons with Br)
        group = "hfo"
        name = name.lower()
        atom_info = name.split(group)[1]
        def_num = "".join([char for char in atom_info if char.isdigit()])
        [no_unsaturated_bonds, no_carbs, no_hydrs, no_fluors] = [int(i) for i in def_num]
        name_PubChem = name.replace("hfo","HFO-")

    elif name.lower().startswith("cfc") or name.lower().startswith("hfc") or name.lower().startswith("hcfc") or name.lower().startswith("hbfc"):
        # chlorofluorocarbons (CFC: C, F, Cl), hydrofluorocarbons (HFC: H, F, C), hydrocarbofluorocarbons (HCFC: H, C, F, Cl), hydrobromochlorocarbons (HBFC: H, Br, F, C)
        grouplist = [char for char in name if char in ["C","F","H","B"]]
        group = "".join(grouplist)
        atom_info = name.split(group)[1]
        def_num_list = [char for char in atom_info if char.isdigit()]
        #use definition of nomen clature, add 90 and calc. number of atoms
        def_num = str(int("".join(def_num_list)) + 90)
        [no_carbs, no_hydrs, no_fluors] = [int(i) for i in def_num]
        # the rest of the non occupied bonds will be filled up with chlorines
        no_chlorines = 2*no_carbs + 2 - no_hydrs - no_fluors
        name_PubChem = name.replace(group, group.upper() + "-")

    elif name.lower().startswith("hcfo"):
        # chlorofluorocarbons (CFC: C, F, Cl), hydrofluorocarbons (HFC: H, F, C), hydrocarbofluorocarbons (HCFC: H, C, F, Cl), hydrobromochlorocarbons (HBFC: H, Br, F, C)
        grouplist = [char for char in name if char in ["C","F","H","B"]]
        group = "".join(grouplist)
        atom_info = name.split(group)[1]
        def_num_list = [char for char in atom_info if char.isdigit()]
        #use definition of nomen clature, add 90 and calc. number of atoms
        def_num = str(int("".join(def_num_list)) + 90)
        [no_unsaturated_bonds, no_carbs, no_hydrs, no_fluors] = [int(i) for i in def_num] 
        # the rest of the non occupied bonds will be filled up with chlorines, but reduced by each unstaturated/ double binding which consumes two of the 'potentially free' spots more than normal binding.
        no_chlorines = 2*no_carbs + 2 - no_unsaturated_bonds*2 - no_hydrs - no_fluors
        name_PubChem = name.replace(group, group.upper() + "-")

    else:
        if isinstance(Formula(name), Formula):
            [formula, name_PubChem] = [name, ""]
        else:
            print("no formula nor classical halocarbon nor in exception list, add formula for substance manually to exception list")
    atnames = ["C","H","Br","Cl", "F"] #Hill notation: C, H, others in alphabetical order
    formula = ""
    for i, el in enumerate([no_carbs, no_hydrs, no_bromines, no_chlorines, no_fluors]): 
        if el==1:
            formula_part = atnames[i]
            formula += formula_part
        elif el>0: 
            formula_part = atnames[i]+str(el)
            formula += formula_part
    return formula, name_PubChem

def search_for_compound_in_PubChem_Lib(compound_name):
    """ return smiles code from a automated web search (PubChem) in chrome, no error handling yet"""
    import selenium
    from selenium import webdriver
    # Using Chrome to access web
    from webdriver_manager.chrome import ChromeDriverManager
    #driver = webdriver.Chrome(r"C:\Users\kaho\.wdm\drivers\chromedriver\win32\102.0.5005.61")
    # Open the website
    #driver.get('https://pubchem.ncbi.nlm.nih.gov/')
    #set download directory path
    from selenium.webdriver.chrome.options import Options
    from selenium.common.exceptions import NoSuchElementException
    import requests
    driver = webdriver.Chrome(ChromeDriverManager().install())

    options = Options()
    options.add_argument("--headless") #without showing the browser
    options.add_experimental_option("prefs", {
    "download.default_directory": r"C://Users//kaho//Downloads//Test",
    "download.prompt_for_download": False,
    "download.directory_upgrade": True,
    "safebrowsing.enabled": True
    })
    #capabilities['chromeOptions']['prefs']['download.prompt_for_download'] = False
    #capabilities['chromeOptions']['prefs']['download.default_directory'] = self.driver_settings['download_folder']
    driver.implicitly_wait(2) #to make sure everything is loaded
    #    search_box = driver.find_element_by_css_selector("[id^='search_']")
    #    search_box.send_keys('\"HFC-152a\"')
    #    search_box.submit()
    #compound_name = 'HFC-152a'
    search_str = 'https://pubchem.ncbi.nlm.nih.gov/#query=%22'+ compound_name +'%22&tab=compound'
    driver.get(search_str)
    #WebDriverWait(driver, 10)
    try:
        search_box = driver.find_element_by_id("Download")
        flag = True
    except NoSuchElementException:
        flag = False
    if flag:
        search_box.click()
        result_box = driver.find_element_by_id("collection-results-container")
        dl_json_buttom = result_box.find_element_by_css_selector("[data-label^='Download; JSON - Save']")
        href_summary = dl_json_buttom.get_attribute('href')

        response = requests.get(href_summary, timeout=15)
        data = response.text

        smiles_str_h = data.split("\"isosmiles\": ")[1].split(",")[0]
        smiles = str(smiles_str_h.replace("\'","").replace("\"",""))
        print(smiles)
        return smiles
    else:
        return None

#Manual Tools
def files2database(list_of_abs_paths_files):
    """ converts all files of filelist to one mgf file, which have the confirmed substance name within the filename (between the characters _ and .), e. g. xxxx_HFC134a.txt.
        It converts the name into a formula using conventions for cfc, hfc, hcfc, hbfc and halons and if not interprets the name as a formula, or searches for an entry in the exception list. 
        The function should only be used with manually control, or if it you can be sure, that conditions are met. It is not meant for beeing called in an automized routine (no error handling yet).
        Maybe control by user IO will be added at a later stage. 
    """
    except_dict = {"HFC4310mee": "C5H2F10", "hexachlorobutadiene": "C4Cl6", "TCHFB": "C4Cl4F6", "PFCc318": "C4F8"}
    alpinac_data_base_list = []
    for path in list_of_abs_paths_files:
        # load
        # searches stuff between _ and . starting with a letter and followed by letters and digits, where (?) is non-greedy: shortest str will be shown
        pattern = "\_[a-zA-Z][a-zA-Z0-9]*?\."
        filename = Path(path).name
        #filename = '190424.0810.tank.3.frag.val_HFC134a.comp.1478.65_0.txt'
        x = re.search(pattern, filename);print(x)
        if not x:
            print("Substance could not extracted from filename: "+filename)
            continue
        else:
            substance = x.group().replace("_","").replace(".","")
            
            #substance = Path(path).name.split("_",".")[-1].split("\n")[0] #does not always work
            instance = AlpinacData.from_alpinac(path)
            # add massbank required metadata (see zimwiki for future plans!)
            instance = AlpinacData.from_alpinac(path)
            instance.set('compound_name', substance)
            formula_prec, name_PubChem = halocarbname2formula(substance, except_dict)
            instance.set('precursor_formula', formula_prec)

            instance.set('precursor_mz', chem_formula(formula_prec).isotope.mass)
            smiles_by_search = search_for_compound_in_PubChem_Lib(name_PubChem)
            if name_PubChem:
                instance.set('smiles', smiles_by_search)
            # append to list for mgf files
            alpinac_data_base_list = alpinac_data_base_list + [instance]
       
    save_as_mgf(alpinac_data_base_list, r"C:/Users/kaho/polybox/ALPINAC\test_kaho\test6.mgf")

 


if __name__ == "__main__":
    #JUST PLAYING AROUND HERE, TODO: DELETE FOR FINAL VERSION
    

        # dl_json_buttom.click() automatic in download folder


# TODO download into specified FOLDER (Empa does not seem to permit it)
# or get from href



#   attr_data_label = result_box.get_attribute("data-label")
#   print(attr_data_label)

#   Download_box = driver.findElement(By.xpath ("//*[contains(text(),'data-label=\"Download; JSON - Save\" ')]"));
#    Download_box = @FindBy(xpath = "//div[@data-label='Download; JSON - Save']")



    plot_manually = False # 12.06.23 examples
    # if plot alpinac example data manuelly
    if plot_manually:
        path_result = r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230401.1747.air.3\230401.1747.air.3.frag.1835s1843s_EI_only_min_4_peaks_per_compound\Compound_0\results_file_mygu.txt"
        path_result = r'C:\\Users\\kaho\\Desktop\\data\\data_Empa\\Campaign202303_interesting_file\\230311.0125.tank.1\\alpinac_results\\230311.0125.tank.1_EI_CI_output_1_mod\\Compound_4\\results_file_mygu.txt'
        alp_spec = [AlpinacData.from_alpinac(path_result)]
        fig, ax = plt.subplots(figsize=(8, 6))
        alp_spec[0].plot_stem(ax, label_peaks = "numeric")
        # set title
        ax.set_title(('1835s_to_1843s, cmp_0'))


    def plottest(ax, figsize=(10, 8), **kwargs):
        if not isinstance(figsize, tuple) or len(figsize) != 2 or \
        not all(isinstance(v, (int, float)) and v > 0 for v in figsize):
            raise ValueError("Invalid figsize argument. Please provide a tuple of positive floats.")
    # rest o

    def plottest2(ax, figsize=(10, 8)):
        if not isinstance(figsize, tuple) or len(figsize) != 2 or \
        not all(isinstance(v, (int, float)) and v > 0 for v in figsize):
            # If figsize is invalid, set it to a default value
            figsize = (10, 8)
            print("Invalid figsize argument. Setting to default value (10, 8).")
        # rest of the code



    
    # TODO write a init procedure, that does the processing,
    # creates an hash or use filesize and if this did not change, dont do preprocessing again but loads it directly
    # export in which file format?
    path_data = "C:/Users/kaho/polybox/ALPINAC/test_kaho/test_data"
    file_mgf = os.path.join(path_data,
                            "GNPS-NIH-NATURALPRODUCTSLIBRARY.mgf")
    # spectrums is a Python list of matchms-type Spectrum objects
    spectrums = list(load_from_mgf(file_mgf))
    
    # process meta data, normalize intensities, reject intensities below a limit and select a certain range
    spectrums = [metadata_processing(s) for s in spectrums]
    spectrums = [peak_processing(s) for s in spectrums]
    
    file_export_json = os.path.join(
        path_data, "GNPS-NIH- NATURALPRODUCTSLIBRARY_processed.json")
    save_as_json(spectrums, file_export_json)

    # pickle should be faster
    file_export_pickle = os.path.join(
        path_data, "GNPS-NIH-NATURALPRODUCTSLIBRARY_processed.pickle")
    pickle.dump(spectrums,
                open(file_export_pickle, "wb"))
    
    
    # plot spectrum with matchms lib:
        # Top spectrum.
    fig, ax = plt.subplots(figsize=(8, 6))
    # plot_spectrum(spectrums[0], mirror_intensity=False, ax=axes, peak_color="C0")
    #y_max = ax.get_ylim()[1]

    plot_mass_spectra(ax, spectrums[0], 
                          spectrums[1],
                          no_plt_legend = 1)
    plt.show()

    #instance = AlpinacData(np.array([1.,2.,3.]), np.array([1.,2.,3.]))
    #instance.uncertainties = np.array([1,2,3])
    path = r"C:\Users\kaho\polybox\ALPINAC\Alpinac_packing\halosearch\data\nontarget_screening\formulas\200602.1550.std.7.frag.train_C2H6.comp.1275.56_0.txt"
    path = r"C:\Users\kaho\polybox\ALPINAC\alpinac\data\nontarget_screening\formulas\190424.0810.tank.3.frag.val_HFC236cb.comp.1556.04_0.txt"
    path = str(Path('C:/Users/kaho/Desktop/data/data_Empa/Campaign202303/230401.1747.air.3/230401.1747.air.3.frag.1647s1661s_EI_only_min_4_peaks_per_compound/Compound_0/results_file_mygu.txt'))
    #path = r'C:/Users/kaho/polybox/ALPINAC/alpinac/data/nontarget_screening/fragments/190424.0810.tank.3.frag.val_HFC236cb.txt'
    instance = AlpinacData.from_alpinac(path)
    #bla = [float(val) for val in instance.metadata['compound_index'][1:-1].split(',')]
    
    #import types
    #a = types.SimpleNamespace()
    plt.rcParams['figure.figsize'] = (8,6)
    # Get type of plt.rcParams['figure.figsize']
    type(plt.rcParams['figure.figsize'])
    # Set the plt.rcParams['figure.figsize'] to a new value
    plt.rcParams['figure.figsize'] = (10, 8)
    # Get type of plt.rcParams['figure.figsize']
    type(plt.rcParams['figure.figsize'])
    

    a = AlpinacData.from_alpinac(path)
    #fig, ax = plt.subplots(1, 1, dpi = 100)
    fig, ax = plt.subplots(figsize=(8, 6))
    #plt.gcf().set_size_inches(3, 2)
    #ax.set_size_inches(8, 6) # change the size of the figure here
    #plottest(ax, plt.rcParams['figure.figsize'])
    plottest2(ax, plt.rcParams['figure.figsize'])

    #plt.plot(a.peaks.mz, a.peaks.intensities, 'o', color = 'C0')
    #stemlines, markerlines, _ = plt.stem(a.peaks.mz, a.peaks.intensities)
    a.plot_stem(ax)

    b = AlpinacData.from_alpinac_unidentified(path)
    #fig, ax = plt.subplots(1, 1, dpi = 100)
    fig, ax = plt.subplots(figsize=(8, 6))
    #plt.gcf().set_size_inches(3, 2)
    #ax.set_size_inches(8, 6) # change the size of the figure here
    #plottest(ax, plt.rcParams['figure.figsize'])

    #plt.plot(a.peaks.mz, a.peaks.intensities, 'o', color = 'C0')
    #stemlines, markerlines, _ = plt.stem(a.peaks.mz, a.peaks.intensities)
    b.plot_stem(ax)


     

    fig, ax = plt.subplots(figsize=(8, 6))
    # plot_spectrum(spectrums[0], mirror_intensity=False, ax=axes, peak_color="C0")
    #y_max = ax.get_ylim()[1]

    plot_mass_spectra(ax, a, 
                          a,
                          no_plt_legend = 0)
    plt.show()

    #b = AlpinacData.from_alpinac_unidentified(path)
    fig, ax = plt.subplots(1, 1)
    
    #b.plot(ax)
    plot_mass_spectra(ax, a, a)


    #type(a) == ms.Spectrum.Spectrum
    #a = spectrums[0]
    best_matches, scores, no_fitting_peaks = score_matching(a, ref_spectra = [a],
                                    tolerance_level = 0.005, min_match_peaks = 5,
                                    no_best_matches = 10)
    # Top spectrum.
    fig, axes = plt.subplots(figsize=(1, 1))
    ax = axes.Axes()
    best_matches[0].plot(ax)
    type(best_matches[0])
    
    database = list(ms.importing.load_from_mgf(r"C:\Users\kaho\polybox\ALPINAC\Alpinac_packing\halosearch\alpinac_gui\data\test5.mgf"))
    [spectra, score, matching_peaks] = score_matching(a, database, tolerance_level=0.01)
    fig, axes = plt.subplots(1, 1)
    ax = axes
    plot_mass_spectra(ax, a, spectra[0],1)
    plt.show()
    
    plot_mass_spectra(ax, a, spectra[1])
    AlpinacData.peaks
    type(a)
    a.metadata
    #a.precursor_formula
    new_var = a.__getattribute__
    new_var
    a.peaks.mz[0]
    #fig, axes = plt.subplots(1, 1)
    #ax = axes
    a.plot(ax)
#    a.keys()
    a.__dict__
    
    name = 'HFC236cb'
    name = 'CFC11bb'
    #following info at https://cdiac.ess-dive.lbl.gov/pns/cfcinfo.html

    name="Halon 1211"       
    halocarbname2formula("Halon 1211")
    
    from typing import List
    import pyteomics.mgf as py_mgf
    from matchms.Spectrum import Spectrum
    a.ionisation_mode
    
    #from molmass import Formula as chem_formula
    #calculate isotopic spectrum
    f = chem_formula(a.metadata.get("est_precursor_formula").split(",")[0])
    f = chem_formula('H2C3F6')
    f = chem_formula('H')
    f = chem_formula('C')
    #f = chem_formula('HFC236cb')

    f.empirical
    f.mass #average mass
    f.isotope.mass
    print(f.atoms)
    print(f.composition())
    s = f.spectrum()
    print(f.spectrum())
    
    save_as_mgf(a, r"C:/Users/kaho/polybox/ALPINAC\test_kaho\test.mfg")

    #C3H2F6 ; Molar Mass, 152.03
   

    path = r"C:\Users\kaho\polybox\ALPINAC\Alpinac_packing\halosearch\data\nontarget_screening\formulas\200602.1550.std.7.frag.train_C2H6.comp.1275.56_0.txt"
    
    
    # convert Myriams file into mfg
    import glob
    #print(files)
    
    # convert identified spectra from Myriam into database
    list_of_abs_paths_files = glob.glob(r"C:\Users\kaho\polybox\ALPINAC\alpinac\data\nontarget_screening\formulas\[2]*")
    files2database(list_of_abs_paths_files)
    
    
    alpinac_data_base_list = []
    

# TODO: sth with the graphical representation is for sure not correct
# plot function does not always work
# why has Myriams file such a weird multipeak shape? Do I use the correct intensity?
        
        
        # TODO decide what to add later: T, gaschromatography column and flux determine RT, so they are important
        # ['pepmass', 'charge', 'mslevel', 
        # 'source_instrument', 'seq', 'ionmode',
        # 'organism', 'smiles', 'inchi',
        # 'inchiaux', 'pubmed', 'libraryquality',
        # 'scans', 'file_name', 'compound_name',
        # 'principal_investigator', 'data_collector',
        # 'submit_user', 'spectrum_id', 'precursor_mz',
        # 'adduct', 'inchikey']

    search_for_compound_in_PubChem_Lib("Halon-123a")