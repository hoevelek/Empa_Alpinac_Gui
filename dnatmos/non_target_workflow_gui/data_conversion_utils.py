import numpy as np
from data_analysis_utils import get_chemical_possible_combinations
from molmass import Formula as chem_formula
from rdkit import Chem
from alpinac.io_tools import get_data_from_jdx_file
from jcamp import jcamp_read #to read NIST spectra



def unit_mass_to_pseudo_high_res_spec(mass:np.array, int:np.array):
    """
    Converts unit mass resolution data to pseudo high resolution data by
    comparison to most likely formulas.
    :param mass: mass data
    :param int: intensity data
    :return: high_res_mass, int
    """


# if main
if __name__ == "__main__":
    # Jcamp is not a mandatory import for us
    filename = r"C:\Users\kaho\Desktop\data\data_Empa\tests\75-09-2-Mass.jdx"
    with open(filename, 'r') as f:
        dict_NIST = jcamp_read(f)

    #frag_data_batch = {0:[]}


    NIST_formula = chem_formula(''.join(dict_NIST['molform'].split()))
    NIST_mass = NIST_formula.isotope.mass
    # get all atoms in the formula
    no_atoms = NIST_formula.atoms
    composition = NIST_formula.composition()
    elements = [str(atom_comp) for atom_comp in composition]
    composition['C'].count



    mz = dict_NIST['x']
    int = dict_NIST['y']

    # get possible candidates for highest mass
    # get highest mz and int
    mz_max = mz[np.argmax(int)]
    int_max = np.max(int)

    candidates_gs, candidates_rad = get_chemical_possible_combinations(elements = elements, target_masses = [mz_max], radicals_allowed=True)
    candidates = list(set(candidates_gs + candidates_rad))
    print(chem_formula('ClCH2').isotope.mass)
    # delete candidates that are not possible: e.g. if possible candidate is <= NIST_formula:
    possible_candidates = [cand for cand in candidates if all([cand.composition()[elem].count < NIST_formula.composition()[elem].count for elem in [str(atom_comp) for atom_comp in cand.composition()]])]
    cand_mayor_iso = possible_candidates[0]
    # get the full isoptopic spectrum
    isotope_unit_mass = [key for key in cand_mayor_iso.spectrum().keys()]
    [val.mz for val in cand_mayor_iso.spectrum().values()]
    rel_ints = [val.fraction for val in cand_mayor_iso.spectrum().values()]





    # frag_metadata = NISTMetadata(
    #     chem_name = dict_target['title'],
    #     target_formula = target_formula,
    #     mol_mass_exact = formula_to_mass(target_formula),
    #     CAS = dict_target['cas registry no']  if 'cas registry no' in dict_target else 'no_CAS_' + dict_target['title']
    # )
    # data_NIST, metadata_NIST = get_data_from_jdx_file()
    # data_NIST = data_NIST[0]
