from email import charset
import io
import logging
from pathlib import Path
import tempfile
from fpdf import FPDF
from PIL import Image
import pandas as pd
import numpy as np
from alpinac_gui.matchms_funcs import formula_from_namestring, formula_to_latex, halocarbname2formula, plot_mass_spectra, load_from_mgf, metadata_processing, peak_processing, AlpinacData, StickSpectra, halocarbname2formula
from pyvalem.formula import Formula
from molmass import Formula as chem_formula
import pubchempy as pcp
from matplotlib import pyplot as plt
from matplotlib import axes
from db_reader import AlpinacFragmentsReader, JdxReader, AbstractReader
#from JdxReader import read_file
from matchms.Spectrum import Spectrum
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import chardet
# import JdxReader
import jcamp
from mendeleev import element
from data_analysis_utils import combinations_sum, get_chemical_possible_combinations_for_elem_list
import urllib.request
import os
from data_analysis_utils import get_chemical_possible_combinations
import cirpy
from matchms.exporting import save_as_mgf
from matchms import Spectrum

DBNIST_WEBPATH = "https://webbook.nist.gov/cgi/cbook.cgi?JCAMP="
download_folder_NIST_pseudo_high_res = r"C:\Users\kaho\Desktop\data\data_Empa\NIST_data\pseudo_high_res_NIST"
download_folder_NIST = r"C:\Users\kaho\Desktop\data\data_Empa\NIST_data\unit_mass_NIST"
path_to_possible_fragments = r"C:\Users\kaho\Desktop\data\data_Empa\combinatorial_data"

# test if folders exist ans create if not
if not os.path.exists(download_folder_NIST_pseudo_high_res):
    os.makedirs(download_folder_NIST_pseudo_high_res)
if not os.path.exists(download_folder_NIST):
    os.makedirs(download_folder_NIST)
if not os.path.exists(path_to_possible_fragments):
    os.makedirs(path_to_possible_fragments)

def get_ms_data_from_NIST(cas, index: int = 0, download_folder: str = download_folder_NIST, dbNISTpath: str = DBNIST_WEBPATH ):
    """
    This function downloads the NIST spectra for a given cas number from the NIST database.
    cas: cas number of the compound
    index: index of the spectrum (default = 0)
    download_folder: path to the folder where the spectra are downloaded to
    dbNISTpath: path to the NIST database
    return: NIST spectra
    """

    # download the NIST data from the NIST database (online)from the path https://webbook.nist.gov/cgi/cbook.cgi?ID=C106467&Units=SI&Mask=200#Mass-Spec
    downloadNISTspecpath = str(dbNISTpath + "C" + str(cas) + "&Index=" + str(index) + "&Type=Mass")
    # Create the directory if it does not exist
    if not os.path.exists(download_folder):
        os.makedirs(download_folder)
    # Extract the filename from the URL
    filename = downloadNISTspecpath.split('/')[-1].split("JCAMP=C")[1].split("&Type")[0].replace("&Index="+str(index),".jdx")
    # Check if the file already exists in the directory
    if not os.path.isfile(os.path.join(download_folder, filename)):
        # Download the file and save it to the specified directory
        urllib.request.urlretrieve(downloadNISTspecpath , os.path.join(download_folder, filename))
        print('File downloaded successfully.')
    else:
        print('File already exists in the directory.')
     # load the file with jcamp
    jdx_reader = JdxReader(download_folder)
    #jdx_reader.read_file(os.path.join(download_folder, filename))
    spectra = jdx_reader.read_spectra(str(Path(filename).stem))
    return spectra

def isotopologues_spectrum(most_abundant_isotopologue_formula):
    """
    This function calculates the isotope pattern for a given formula.
    most_abundant_isotopologue_formula: formula of the most abundant isotopologue
    return: mass_isotope_formula: dictionary with the masses as keys and the formulas as values
    """
    # calculate the isotope pattern for the constitutents of the fragment
    atoms = [at for at in chem_formula(most_abundant_isotopologue_formula).composition().asdict().keys()]
    isotopes_all = []
    dict_elem_of_isotopes_all = {}
    #dict_dm_of_isotopes_to_most_abundant = {}
    # calcualte differences beween isotopes to be able to assign them to the correct element
    for at in atoms:
        atom_isotopic_spectra = chem_formula(at).spectrum().asdict()
        isotope_mass_atom = [mass for mass in atom_isotopic_spectra.keys()]
        isotopes_all.extend(isotope_mass_atom)
    # isotope_mass_atom_exact = [atom_isotopic_spectra[mass][0] for mass in chem_formula(at).spectrum().asdict().keys()]
        #dm_isotopes = abs(chem_formula(at).monoisotopic_mass - np.array(isotope_mass_atom_exact))

        # add the element to the dictionary
        dict_elem_of_isotopes_all.update({mass: "[" + str(mass) + str(at) + '], ' + str(-chem_formula(at).monoisotopic_mass + atom_isotopic_spectra[mass][0]) for mass in isotope_mass_atom})
        # calculate the differences between all the isotopes amongst each other
        # delete zero
        #dm_isotopes = np.delete(dm_isotopes, np.where(dm_isotopes == 0))
    # calculate the isotope pattern for the fragment
    mass_isotope_pattern = chem_formula(most_abundant_isotopologue_formula).spectrum().asdict()
    mass_isotope_formula = {}
    # get the largest isotope peak
    most_abundant_isotopologue_mass = chem_formula(most_abundant_isotopologue_formula).spectrum().peak.mass
    # eleimiate all with intensity < 0.001
    mass_isotope_pattern = {key: value for key, value in mass_isotope_pattern.items() if value[1] > 0.0001}
    # get the masses of the isotope spectrum
    mass_isotope_pattern.keys()
    # calculate the possible fragments
    for mass in mass_isotope_pattern.keys():
        formulas = []
        # only chemical possible candidates
        #elements_for_chem_formula = [dict_elem_of_isotopes_all[iso] for iso in isotopes_all]
        dict_formula_of_dm = {iso:dict_elem_of_isotopes_all[iso]  for iso in isotopes_all}
        # get the dm of the isotopes to the most abundant isotope # mal 1e6
        dm_isotopes = [int(float(dict_formula_of_dm[iso].split(',')[1])*1e4) for iso in isotopes_all]
        # get the corresponding isotopes
        symbols = [dict_formula_of_dm[iso].split(',')[0] for iso in isotopes_all]
        dict_symbol_of_mass = {iso:symbols[i] for i, iso in enumerate(dm_isotopes)}
        # remove dm_isotopes that are zero
        candidates_masses = [dmi for dmi in dm_isotopes if dmi != 0]
        # convert to list if not already
        if type(candidates_masses) != list:
            candidates_masses = [candidates_masses]
        # get the delta masses of the isotopes to the most abundant isotope
        target_dm = int((mass_isotope_pattern[mass][0] - most_abundant_isotopologue_mass)*1e4)
        # if next to zero, candidates = "most_abundant"
        if abs(target_dm) < 1:
            candidates = [0]
            formula = most_abundant_isotopologue_formula
            formulas.append(formula)
            print(formulas)
        else:  
            candidates = combinations_sum(target = target_dm, candidates = candidates_masses)
            # handling rounding error/ this is not really nice, because the real mass is not the sum of the atom masses!
            nreps = 1
            while len(candidates) == 0 and nreps <= 2:
                # if no candidates are found, increase the target mass by 1
                candidates = combinations_sum(target = target_dm + 1*nreps, candidates = candidates_masses)
                if len(candidates) == 0:
                    candidates = combinations_sum(target = target_dm - 1*nreps, candidates = candidates_masses)
                nreps += 1

        ##candidates_gs, candidates_rad = get_chemical_possible_combinations(elements = elements_for_chem_formula, target_masses = mass, radicals_allowed=True)   
        #candidates = list(set(candidates_gs + candidates_rad))
        
            for cand in candidates:
                if cand == [0]: # formula is most abundant isotope
                    formula = most_abundant_isotopologue_formula
                    formulas.append(formula)
                    print(formulas)
                else: # formula has one or more atoms substituted by isotopes
                    # get the formula of the candidate
                    formula = [dict_symbol_of_mass[dm] for dm in cand]
                    # add all elements together
                    formula_str = ''.join(formula)
                    # remove numbers and [] from the formula
                    formula_main_isotopes_str  = ''.join([i for i in formula_str if not i.isdigit()]).replace('[', '').replace(']', '')
                    #for formula_atom in formula:
                        #formula_atom_only_charaters = ''.join([i for i in formula_atom if not i.isdigit()]) 
                    formula_isotopologue = (chem_formula(most_abundant_isotopologue_formula) - chem_formula(formula_main_isotopes_str) + chem_formula(formula_str)).formula
                    formulas.append(formula_isotopologue)
                    # for each element in the formula, get the corresponding symbol
                    #for i, elem in enumerate(formula):
                    #    formula[i] = dict_symbol_of_mass[elem]
                    # get the formula of the candidate
                    #formulas.append(''.join(formula))
                    print(formulas) 
        if len(formulas) == 0:
            mass_isotope_formula[mass] = ["?"]             
        else:
            mass_isotope_formula[mass] = formulas
    return mass_isotope_formula, mass_isotope_pattern

def get_NIST_EI_mass_spectra_by_identifier(identifier, identifier_type = 'cas', download_folder_local = None, path_to_NIST_database = DBNIST_WEBPATH, plot = False, calculate_high_res: bool = False):
    """
    This function downloads the NIST spectra for a given identifier (cas, inchi, smiles, name, formula, cid, canonical_smiles, smiles) from the NIST database.
    if calculate_high_res is True, the high resolution spectra are estimated from the unit mass spectra by fitting the isotopic spectra of assigned main isotope to it.
    identifier in formate: cas, inchi, smiles, name, formula, cid, canonical_smiles, smiles
    identifier_type: cas, inchi, smiles, name, formula, cid, canonical_smiles, smiles
    download_folder_local: path to local folder where the NIST spectra are downloaded to
    path_to_NIST_database: path to the NIST database
    plot: if True, the NIST spectra are plotted
    calculate_high_res: if True, the high resolution spectra are estimated from the unit mass spectra by fitting the isotopic spectra of assigned main isotope to it.  
    """
    cas_no = []
    # generate NIST spectrum
    if identifier_type == 'inchi':
        cas_no = cirpy.resolve(identifier, 'cas')
        inchi = identifier
    elif identifier_type == 'canonical_smiles':
        cmp=pcp.get_compounds(identifier, 'smiles') #of cid, name, smiles, sdf, inchi, inchikey or formula.
        # get formula from smiles
        sum_formula = [chem_formula(cmp_i).formula for cmp_i in cmp]
        cas_no = [cirpy.resolve(cmp_i.inchi, 'cas') for cmp_i in cmp]
        inchi = [cmp_i.inchi for cmp_i in cmp]
    elif identifier_type == 'smiles':
        cas_no = cirpy.resolve(identifier, 'cas')
        inchi = cirpy.resolve(identifier, 'inchi')
    elif identifier_type == 'name':
        cas_no = []
        try: 
            cas_no = cirpy.resolve(identifier, 'cas')
            inchi = cirpy.resolve(identifier, 'inchi')
        except:
            cmp = pcp.get_compounds(identifier, 'name')
            inchi = [cmp[i].inchi for i in range(len(cmp))]
            cas_no = [cirpy.resolve(inchi, 'cas') for inchi in inchi]
    elif identifier_type == 'cid':
        cmp=pcp.get_compounds(identifier, 'cid')
        cas_no = [cirpy.resolve(cmp_i.inchi, 'cas') for cmp_i in cmp]
        inchi = [cmp_i.inchi for cmp_i in cmp]
    elif identifier_type == 'formula':
        cmp=pcp.get_compounds(identifier, 'formula')
        # print number and name of found names
        names = [cmp_i.iupac_name for cmp_i in cmp]
        inchi = [cmp_i.inchi for cmp_i in cmp]
        print("name found for formula: "+ identifier + str(names))
        print("name found: " + str(len([cmp_i.iupac_name for cmp_i in cmp if cmp_i.iupac_name != None])))
        print('inchi found for formula: ' + identifier + str(inchi))
        print("inchi found: " + str(len([cmp_i.inchi for cmp_i in cmp if cmp_i.inchi != None])))
        #remove None inchis and names
        names = [name for name in names if name is not None]
        inchi = [inchi_i for inchi_i in inchi if inchi_i is not None]
        cas_no = [cirpy.resolve(inchi_i, 'cas') for inchi_i in inchi]
        print("cas_no found for " + str(len([cas for cas in cas_no if cas != None])) + ' from total of ' + str(len(cas_no))+ 'inchis')
    elif identifier_type == 'cas':
        cas_no = [identifier]

        #cas_no_names = [cirpy.resolve(name, 'cas') for name in names] # maybe inchi more suited, for sure faster
        #len(np.unique([str(cs) for cs in cas_no_names]))
        #print("cas_no found for " + str(len([cas for cas in cas_no_names if cas != None])) + ' from total of ' + str(len(cas_no_names))+ 'names')       
    else:
        # raise error: identifier type not known
        print("identifier type not known")


    if identifier_type in ['inchi', 'canonical_smiles',  'name', 'cid', 'formula']:
        #inchi = [cmp_i.inchi for cmp_i in cmp]
        pass


    if identifier_type in ['smiles', 'cas']:
        # determine inchi for all cas_no
        #inchi = []
        #for list_elem in cas_no:
        #    if isinstance(list_elem, list):
        #        inchi.append([cirpy.resolve(cas, 'inchi') for cas in list_elem])
        #    else:
        #        inchi.append(cirpy.resolve(list_elem, 'inchi'))
        try:
            inchi = [cirpy.resolve(cas, 'inchi') for cas in cas_no]
        except:
            inchi = [None]
        #cas = cas_no
        # resolve cas_no to inchi
        #cas_elems = [cas for cas in cas_no if cas != None]
        #cas_test = [cas for cas in cas_elems[0]][0] # does not matter which cas_no is taken, because same formula
        #inchi = cirpy.resolve(cas_test, 'inchi')
        if inchi[0] == None:
            print("no inchi found for cas_no " + cas_no[0])
            cmp = [None]
            #raise error
        else:
            cmp= pcp.get_compounds(inchi, 'inchi')
        
        # get atoms, not matter which cmp is taken, because same formula
        

    inchi = [inchi[i] for i, cas in enumerate(cas_no) if cas != None]
    cas_no = [cas for cas in cas_no if cas != None]
    cmp = [cmp[i] for i, cas in enumerate(cas_no) if cas != None]
    if cmp[0] != None:
        sum_formula = cmp[0].molecular_formula #same for all
    else:
        sum_formula = None


    #TODO DICT INCHI CAS_NO
    inchi_cas_dict = dict(zip(inchi, cas_no))
    cas_inchi_dict = dict(zip(cas_no, inchi))


    #cmp_ind = 0
    #index = str(0)
    #cmp=pcp.get_compounds('benzene', 'name')
    #download_folder = r"C:\Users\kaho\Desktop\data\data_Empa\NIST_data\pseudo_high_res_NIST"
    #cid  = cmp[cmp_ind].cid
    #cas_no  = cirpy.resolve(cmp[cmp_ind].iupac_name, 'cas')
    #cas_no = ['100-42-5']

    if download_folder_local and not os.path.exists(download_folder_local):
        os.makedirs(download_folder_local)

    spectras = []
    superpos_high_res_formulas_dict = {}
    for cas_ind in range(len(cas_no)):
        spectra = get_ms_data_from_NIST(cas_no[cas_ind], index = 0, download_folder = download_folder_local)
        if len(spectra.mz) > 0:
            print(spectra.mz)
            print('spectrum not found in NIST database')
            NIST_mass = spectra.mz
            NIST_int = spectra.intensities/9999
            # to pandas
            NIST_data = pd.DataFrame({'mass': NIST_mass, 'int': NIST_int})
            NIST_stick_spectrum = StickSpectra(mz = NIST_mass, intensities = NIST_int, metadata = {"precursor_mz": -1})
            # set metadata
            NIST_stick_spectrum.set("elements", np.unique(cmp[cas_ind].elements))
            NIST_stick_spectrum.set("formula", cmp[cas_ind].molecular_formula)
            NIST_stick_spectrum.set("inchi", cmp[cas_ind].inchi)
            NIST_stick_spectrum.set("iupac_name", cmp[cas_ind].iupac_name)
            NIST_stick_spectrum.set("cas_no", cas_no[cas_ind])
            NIST_stick_spectrum.set("smiles_canonical", cmp[cas_ind].canonical_smiles)
            isomeric_smiles = cirpy.resolve(cmp[cas_ind].inchi, 'smiles')   
            NIST_stick_spectrum.set("smiles_isomeric", cmp[cas_ind].isomeric_smiles)    

            if plot == True:
                fig, ax = plt.subplots()
                NIST_stick_spectrum.plot_stem(ax, label_peaks = "numeric", label = ["NIST", ""])
                # set title to iupac name and cas no
                plt.title(cmp[cas_ind].iupac_name + ", cas:" + str(cas_no[cas_ind]))   
            # title is the name of the compound
                #plt.title(cmp[cas_ind].iupac_name)
                plt.show()
                #save plot
                if download_folder_local is not None:
                    path_fig = Path(download_folder_local)/Path('plots') /Path(cmp[cas_ind].iupac_name + ".png")
                    # check if folder exists and create if not
                    if not os.path.exists(path_fig.parent):
                        os.makedirs(path_fig.parent)
                    fig.savefig(path_fig)
                spectras.append(NIST_stick_spectrum)
        else:
            # go to next spectrum
            next
            # check if download_folder path exists on computer:

        pseudo_high_res = []
        pseudo_unit_mass = []
        all_isotopologue_specs_all = [] 
        cmp_all = []    
        if (calculate_high_res == True) & (len(spectra.mz) > 0):
            superpos_spectrum, superpos_spectrum_high_res, all_isotopologue_specs, _, _, superpos_high_res_formulas_dict = calculate_high_res_spec(NIST_stick_spectrum)
            pseudo_high_res.append(superpos_spectrum_high_res)
            pseudo_unit_mass.append(superpos_spectrum)
            all_isotopologue_specs_all.append(all_isotopologue_specs)
            cmp_all.append(cmp[cas_ind])    
            # plot the high resolution spectrum
            if plot == True:
                fig, ax = plt.subplots()
                superpos_spectrum_high_res.plot_stem(ax, label_peaks = "numeric", label = ["NIST", ""])
                # set title to iupac name and cas no
                plt.title(cmp[cas_ind].iupac_name + ", cas:" + str(cas_no[cas_ind]))

    return spectras, cmp, inchi_cas_dict, cas_inchi_dict, pseudo_high_res, pseudo_unit_mass, all_isotopologue_specs_all, superpos_high_res_formulas_dict

def calculate_high_res_spec(NIST_spectrum):
    """
    This function calculates the high resolution spectra from the unit mass spectra by fitting the isotopic spectra of assigned main isotope to it.
    NIST_spectrum: NIST spectrum
    return: superpos_spectrum: superposition of all isotopologue spectra
    return: superpos_spectrum_high_res: superposition of all isotopologue spectra in high resolution
    return: all_isotopologue_specs: all isotopologue spectra
    return: cmp: compound  
    return: superpos_high_res_formulas: formulas of the high resolution spectrum

    
    """
    # get elements from spectra metadata
    elements = NIST_spectrum.metadata["elements"]
    # get the formula from the spectra metadata
    formula = NIST_spectrum.metadata["formula"]
    # get the composition of the formula
    composition_dict = chem_formula(formula).composition().asdict()

    # get target masses from the spectra
    NIST_data = pd.DataFrame({'mass': NIST_spectrum.mz, 'int': NIST_spectrum.intensities})
    #target_masses = spectra.mz
    _, _, formulas, masses = get_chemical_possible_combinations_for_elem_list(elements, target_masses=NIST_spectrum.mz, radicals_allowed=True, path_to_files=path_to_possible_fragments, recalculate=False)
    #data_high_res_mass_est = pd.read_csv(high_res_cand, sep="\t", header=None, names = ["fragments", "HR_mass"])
    # order masses, formulas by masses
    masses, formulas = zip(*sorted(zip(masses, formulas)))
    # write into pandas dataframe
    data_high_res_mass_est_h = pd.DataFrame({'fragments': formulas, 'HR_mass': masses})

    for ind, formula in enumerate(data_high_res_mass_est_h['fragments']):
        fragment_comp_dict = chem_formula(formula).composition().asdict()
        # remove rows if formula is not in the composition
        for key in fragment_comp_dict.keys():
            # calculate the isotope pattern for the fragment:
            # check if fragment is valid, meaning containing maximal atoms present in the molecule
            if fragment_comp_dict[key][0] > composition_dict[key][0]:
                print("fragment not valid for atom " + key + " , as it contains more atoms as present in the molecule")
                # set corresponding row to nan
                data_high_res_mass_est_h['HR_mass'].iloc[ind] = np.nan
                data_high_res_mass_est_h['fragments'].iloc[ind] = np.nan
    
    # remove lines containing nan
    data_high_res_mass_est = data_high_res_mass_est_h.dropna()

    # convert data_high_res_mass_est['HR_mass'] to integers
    data_high_res_mass_est['unit_mass'] = np.round(data_high_res_mass_est['HR_mass']).astype(int)
    # find duplicates (radicals + ions)
    data_high_res_mass_est[data_high_res_mass_est.duplicated(subset=['fragments'], keep=False)]
    # sort NIST_data by 'int'
    NIST_data_sorted = NIST_data.sort_values(by=['int'], ascending=False)

    # go through the NIST_data_sorted and calculate isotope pattern for each mass
    NIST_updated_most_abundant_isotopologues = {}
    NIST_isotopologues = {}

    all_isotopologue_specs = []
    for i in range(len(NIST_data_sorted)):
        # get the corresponding fragment data_high_res_mass_est['unit_mass'] to the NIST_data_sorted['mass']
        subset = data_high_res_mass_est[data_high_res_mass_est['unit_mass'] == NIST_data_sorted['mass'].astype(int).iloc[i]]
        # get the corresponding fragment data_high_res_mass_est['unit_mass'] to the NIST_data_sorted['mass']
        for ind_frag_per_peak in range(len(subset)):
            fragment_estimated = subset['fragments'].iloc[ind_frag_per_peak]
            mass_isotopologue_formulas, mass_isotopologue_pattern = isotopologues_spectrum(fragment_estimated)
            mz = np.array([mass_isotopologue_pattern[key][0] for key in mass_isotopologue_pattern.keys()])
            intensities = np.array([mass_isotopologue_pattern[key][2]/100 for key in mass_isotopologue_pattern.keys()])
            formulas = [mass_isotopologue_formulas[key][0] for key in mass_isotopologue_pattern.keys()]
            isotopologue_spectrum = StickSpectra(mz = mz, intensities = intensities, metadata = {"precursor_mz": -1})
            isotopologue_spectrum.set("formula", formulas)
            isotopologue_spectrum.set("main_isotopologue", fragment_estimated)
            all_isotopologue_specs.append(isotopologue_spectrum)


    #iso_stick_spectra = all_isotopologue_specs
    #fig, ax = plt.subplots()
    #spec_to_plot = AlpinacData(mz =  all_isotopologue_specs[0].peaks.mz, intensities= all_isotopologue_specs[0].peaks.intensities,  metadata = {"precursor_mz": -1})
    #spec_to_plot.set("formula",  all_isotopologue_specs[0].get("formula"))
    #spec_to_plot.formula =  all_isotopologue_specs[0].get("formula")  
    ##spec_to_plot.plot_stem(ax, label_peaks = "numeric")
    #spec_to_plot.plot_stem(ax, label_peaks = "formulas")
    # linear dependencies
    (
        lin_dep_spectra_h,
        lin_independent_spectra_h,
    ) = StickSpectra.linear_dependence_of_spectra(
        # tollevel is multiplied by 2, because theo. peaks p +- tollevel []: [p1][p2] could lead to contribution of same measured peak.
        # to exclude that, in the extreme case where measured peak is exactly in the middle of both theoretical peaks [p1][m][p2] the
        # right side of p1 plus the left side of m must give the effective overlap
        list_spectra= all_isotopologue_specs,
        tolerance_level=0.45*2,
        threshold_score_for_lin_dep=0,
    )

    # TODO determine later as rest
    lin_dep_spectra = lin_dep_spectra_h[0]
    lin_independent_spectra = lin_independent_spectra_h[0]  

    # 1) independent spectra: just scale as good as possible
    for spec in lin_independent_spectra:
        k_abs = spec.scale_fac_rel_to(NIST_spectrum, tollevel=0.49)
        spec.set("k_abs", k_abs)

    # 2) dependent spectra: fit aech spectra block a) relative and b) absolute to the NIST spectrum
    for i, spec_set in enumerate(lin_dep_spectra):
        # fit relative
        k_rel_i = NIST_spectrum.fit_absolute_stick_spectra_contributions(
        ref_spectra=spec_set,
        tollevel=0.249,
        )

        # fit absolute
        # get superposition of all spectra in specset
        superpos_spectrum = spec_set[0].superpose(spec_set[1:], k_rel_i[0], k_rel_i[1:], tollevel_aggregate=0.249)
        #k_abs_set = superpos_spectrum.scale_fac_rel_to(NIST_stick_spectrum, tollevel=0.249)
        k_abs_set = NIST_spectrum.scale_fac_rel_to(superpos_spectrum, tollevel=0.249)

        k_abs = k_abs_set*np.array(k_rel_i)
        print(k_abs_set)
        print(k_abs)
        for j, spec in enumerate(spec_set):
            spec.set("k_rel", k_rel_i[j])
            spec.set("superposition_id", i)
            spec.set("k_abs", k_abs[j])
            print(k_abs[j])

    # get flat list of all isotopologue spectra
    all_isotopologue_specs = []
    for spec_set in lin_dep_spectra:
        all_isotopologue_specs.extend(spec_set) 
    for spec in lin_independent_spectra:
        all_isotopologue_specs.append(spec)

    # get all k_abs from metadtata
    k_abs_all = np.array([spec.metadata["k_abs"] for spec in all_isotopologue_specs]) 



    # superposition of all spectra
    superpos_spectrum = all_isotopologue_specs[0].superpose(all_isotopologue_specs[1:], k_abs_all[0], k_abs_all[1:], tollevel_aggregate=0.49)
    superpos_spectrum_high_res = all_isotopologue_specs[0].superpose(all_isotopologue_specs[1:], k_abs_all[0], k_abs_all[1:], tollevel_aggregate=0)
    # get all formulas
    formulas = [spec.get('formula') for spec in all_isotopologue_specs]
    # get all high res masses
    masses = [spec.peaks.mz for spec in all_isotopologue_specs]
    # flatten and convert to dict
    formulas_dict = dict(zip( [item for sublist in masses for item in sublist], [item for sublist in formulas for item in sublist]))
    # get all the formulas for masses in high res superpos spectrum
    formulas_superpos_spectrum = [formulas_dict[mass] for mass in superpos_spectrum_high_res.peaks.mz]
    # convert to AlpinacData
    #superpos_spectrum_high_res = AlpinacData(mz = superpos_spectrum_high_res.peaks.mz, intensities = superpos_spectrum_high_res.peaks.intensities, metadata = {"precursor_mz": -1}) 
    # set formulas to metadata
    superpos_spectrum_high_res.set("formula", str(np.array(formulas_superpos_spectrum)))

    return superpos_spectrum, superpos_spectrum_high_res, all_isotopologue_specs, k_abs_all, NIST_spectrum, formulas_superpos_spectrum

def load_pseudo_high_res_spectrum(cas, path_save, largest_N_peaks: int = None):
# try to load from download folder
    if os.path.exists(path_save):
        pseudo_high_res = list(load_from_mgf(path_save))[0]
        # get formulas from metadata
        formulas_superposspectrum_str = pseudo_high_res.metadata['fragment_formulas']
        # convert string(array) to list
        #remove bracket at start and end
        formulas_superposspectrum_str = formulas_superposspectrum_str[1:-1]
        formulas_superposspectrum = formulas_superposspectrum_str.replace("'", "").split(', ')
    else: #generate it
        NIST_spectrum, cmp, inchi_cas_dict, cas_inchi_dict, pseudo_high_res, pseudo_unit_mass, all_isotopologue_specs_all, formulas_superposspectrum = get_NIST_EI_mass_spectra_by_identifier(identifier=cas, identifier_type = 'cas', download_folder_local = download_folder_NIST, path_to_NIST_database = DBNIST_WEBPATH, plot = True, calculate_high_res=True)
        # convert all metadata to strings
        pseudo_high_res = pseudo_high_res[0]
        meta_pseudo_HR = {'fragment_formulas': str(formulas_superposspectrum)}
        spec_ms = Spectrum(mz = pseudo_high_res.peaks.mz, intensities = pseudo_high_res.peaks.intensities, metadata = meta_pseudo_HR, metadata_harmonization=False) 
        save_as_mgf(spectrums = spec_ms, filename = path_save)
    
    # set formulas (for plotting, ...)
    pseudo_high_res.set("formula", formulas_superposspectrum)
    # get largest 12 peaks and corresponding formulas
    if largest_N_peaks == None:
        pseudo_high_res_largest_N = pseudo_high_res
        formulas_largest_N = formulas_superposspectrum
    else:
        # get indidcies of largest 12 peaks
        ind_largest_N = np.argsort(pseudo_high_res.peaks.intensities)[-largest_N_peaks:]
        # get mz and intensities of largest 12 peaks
        mz_largest_N = pseudo_high_res.peaks.mz[ind_largest_N]
        int_largest_N = pseudo_high_res.peaks.intensities[ind_largest_N]
        # get formulas of largest 12 peaks
        formulas_largest_N = [formulas_superposspectrum[ind] for ind in ind_largest_N]
        # sort by mz
        mz_largest_N, int_largest_N, formulas_largest_N = zip(*sorted(zip(mz_largest_N, int_largest_N, formulas_largest_N)))
        # convert to AlpinacData
        pseudo_high_res_largest_N = AlpinacData(mz = np.array(mz_largest_N), intensities = np.array(int_largest_N), metadata = {"precursor_mz": -1})
        # set formulas to metadata
        pseudo_high_res_largest_N.set("formula", formulas_largest_N)
    return pseudo_high_res, pseudo_high_res_largest_N, formulas_largest_N

if __name__ == "__main__":

    # save spectra


    #get_NIST_EI_mass_spectra_by_identifier(identifier ='C1=CC=CC=C1', identifier_type = 'canonical_smiles', download_folder_local = download_folder_NIST, path_to_NIST_database = DBNIST_WEBPATH, plot = True)
    #get_NIST_EI_mass_spectra_by_identifier(identifier ='Benzene', identifier_type = 'name', download_folder_local = download_folder_NIST, path_to_NIST_database = DBNIST_WEBPATH, plot = True)
    

    cas = "311-89-7"
    cas = "74-85-1"
    cas = '75-02-5'
    path_save = str(Path(download_folder_NIST_pseudo_high_res) / Path(cas + ".mgf"))

    pseudo_high_res, pseudo_high_res_largest_N, formulas_largest_N = load_pseudo_high_res_spectrum(cas, path_save, largest_N_peaks=12)
    fig, ax = plt.subplots()
    pseudo_high_res_largest_N.plot_stem(ax, label_peaks = "formulas")
    #set title to iupac name and cas no
    plt.title(cas)

    cas = '75-02-5'
    #cas = '692-49-9'
    # get NIST spectrum
    NIST_spectrum, cmp, inchi_cas_dict, cas_inchi_dict, pseudo_high_res, pseudo_unit_mass, all_isotopologue_specs_all, formulas_superposspectrum = get_NIST_EI_mass_spectra_by_identifier(identifier=cas, identifier_type = 'cas', download_folder_local = download_folder_NIST, path_to_NIST_database = DBNIST_WEBPATH, plot = True, calculate_high_res=True)

    # plot
    fig, ax = plt.subplots()
    pseudo_high_res[0].plot_stem(ax, label_peaks = "formulas")
    pseudo_high_res[0].metadata['formula']
    # save as mgf
    path_save = str(Path(download_folder_NIST_pseudo_high_res) / Path(cas + ".mgf"))
    save_as_mgf(spectrums = pseudo_high_res[0], filename = path_save)
    # PLOTTING
    
    spec_index = 0
    # plot superpos_spectrum_high_res
    fig, ax = plt.subplots()
    pseudo_high_res[spec_index].plot_stem(ax, label_peaks = "numeric")

    fig, ax = plt.subplots()
    pseudo_high_res[spec_index].plot_stem(ax, label_peaks = "formulas")

    fig, ax = plt.subplots()
    pseudo_unit_mass[spec_index].plot_stem(ax, label_peaks = "numeric")

    fig, ax = plt.subplots()
    NIST_spectrum[spec_index].plot_stem(ax, label_peaks = "numeric")

    # if fit worked, should be one
    scale_global = NIST_spectrum[spec_index].scale_fac_rel_to(pseudo_unit_mass[spec_index], tollevel=0.49)

    k_abs_all = np.array([spec.metadata["k_abs"] for spec in all_isotopologue_specs_all[spec_index]])   


    iso_stick_spectra_alpinac_scaled = [AlpinacData(spec.peaks.mz, k_abs_all[i]*spec.peaks.intensities, metadata=spec.metadata) for i, spec in enumerate( all_isotopologue_specs_all [spec_index])]
    color_list = ["orange", "green", "blue", "red", "purple", "brown", "pink", "gray", "olive", "cyan", "magenta", "yellow", "black", "lightblue", "lightgreen", "lightgray", "lightpink", "lightyellow", "lightcyan"]
    # yaxis from -1 to 1
    #ax.set_ylim([-1, 1]) 
    fig, ax = plt.subplots() 
    for i in range(len(iso_stick_spectra_alpinac_scaled)): 
        #get formulas from metadata
        # get the formulas from the metadata
        formulas = iso_stick_spectra_alpinac_scaled[i].metadata["formula"]
        main_iso = iso_stick_spectra_alpinac_scaled[i].metadata["main_isotopologue"]
        # set formulas to metadata
        iso_stick_spectra_alpinac_scaled[i].plot_stem(ax, col=color_list[i%len(color_list)], label_peaks = "formulas", label = [main_iso,""])

    # add nist spectrum as crosses
    plt.plot(NIST_spectrum[spec_index].mz, NIST_spectrum[spec_index].intensities, 'x', color = 'black', label = "NIST")
    # add superpos_spectrum_high_res as red horizontal - 
    plt.plot(np.round(pseudo_unit_mass[spec_index].mz).astype(int), pseudo_unit_mass[spec_index].intensities, '_', color = 'red', label = "pseudo high res")
    # set title to iupac name and cas no
    plt.title(cmp[spec_index].iupac_name + ", cas:" + str(NIST_spectrum[spec_index].metadata["cas_no"]))
    plt.show()

