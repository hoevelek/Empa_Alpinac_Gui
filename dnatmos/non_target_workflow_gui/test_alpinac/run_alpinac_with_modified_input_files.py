from alpinac.mode_identification.main_from_file import make_identification_from_file
from molmass import Formula as chem_formula
import pubchempy as pcp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

dir_results = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230401.1747.air.3")
ch2cl2_to_theo = False
check_largest_peaks_only = True
test_multiple_factors = False


# open the file C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230401.1747.air.3\230401.1747.air.3.frag.1647s1661s_EI_only_min_4_peaks_per_compound.txt
file = Path(dir_results / '230401.1747.air.3.frag.1647s1661s_EI_only_min_4_peaks_per_compound.txt')
file = Path(dir_results / '230401.1747.air.3.frag.1592s1600s_EI_only_min_4_peaks_per_compound.txt') #CFC-11

dir_results = Path(r"C:\Users\kaho\Desktop\data\data_Empa\tests")
file = Path(dir_results / Path('230401.1747.air.3.frag.1592s1600s.txt'))


#open file as pandas dataframe
df = pd.read_csv(file, sep='\t', header=0, index_col=0)
df = pd.read_csv(file, header=0, delim_whitespace=True, skipinitialspace=True)





if ch2cl2_to_theo == True:
    ## THIS PART MODIFIES INPUT DATA TO THEORETICAL DATA ##
    # # Search for compounds with the name "ethane"
    results = pcp.get_compounds('CCl', 'formula', as_dataframe=True)
    #results = pcp.get_compounds('HCCl', 'formula', as_dataframe=True)
    # # Get the first result from the search
    compound = results[0]
    # # Get the molecular formula from the compound
    formula = compound.molecular_formula
    print(formula)  # Output: C2H6
    mass_exact = pcp.get_properties('ExactMass', compound.cid)
    print(mass_exact)
    # smiles_can = pcp.get_properties('CanonicalSMILES', compound.cid)

    dict_CH2Cl2 = {}

    dict_CH2Cl2[47] = 'CCl'
    dict_CH2Cl2[48] = 'HCCl'
    dict_CH2Cl2[49] = 'H2CCl'
    dict_CH2Cl2[50] = 'H2[13C]Cl'
    dict_CH2Cl2[51] = 'H2C[37Cl]'
    dict_CH2Cl2[52] = 'H2[13C][37Cl]'
    dict_CH2Cl2[82] = 'CCl2'
    dict_CH2Cl2[83] = 'HCCl2'
    dict_CH2Cl2[84] = 'H2CCl2'
    dict_CH2Cl2[85] = 'H2[13C]Cl2'
    dict_CH2Cl2[86] = 'H2CCl[37Cl]'
    dict_CH2Cl2[87] = 'H2[13C]Cl[37Cl]'
    dict_CH2Cl2[88] = 'H2C[37Cl]2'

    # loop over all dictionary values:
    dict_CH2Cl2_exact_mass = {}

    for elem in dict_CH2Cl2:
        print(elem)
        chem_form = chem_formula(dict_CH2Cl2[elem])
        dict_CH2Cl2_exact_mass[elem] = chem_formula(dict_CH2Cl2[elem]).isotope.mass



    # get the column "mass" as a list
    mass_list = df['mass'].tolist()
    # Convert the values of the dictionary to a list
    dict_values = list(dict_CH2Cl2_exact_mass.values())
    # elements of the mass_list are replaced by dict_values choosing the nearest value
    for i in range(len(mass_list)):
        #which index is the nearest value in dict_values to the value in mass_list[i]
        idx = (np.abs(np.asarray(dict_values) - mass_list[i])).argmin()
        # if minimum distance is smaller than 0.1, replace the value in mass_list with the value in dict_values
        if np.abs(np.asarray(dict_values)[idx] - mass_list[i]) < 0.1:
            mass_list[i] = dict_values[idx]
        else:
            mass_list[i] = np.nan
    # replace the column "mass" in the dataframe with the new mass_list
    df['mass'] = mass_list
    # delete all rows with NaN values in the column "mass"
    df = df.dropna(subset=['mass'])

if test_multiple_factors == True:
# multiply the values in the column "mass_u_cal" with factor
    factors = [1.0, 1.25,1.50,1.75,2.0]
    for factor in factors:
        df["mass_u_ppm"] = df["mass_u_ppm"]*factor


        # set all entries of df["Adduct"] to None
        # replace all entries of df["Adduct"]=="NaN" with "None"
        df["Adduct"] = df["Adduct"].replace(np.nan, 'None', regex=True)
        #df["Adduct"] = ["None"]*len(df)
        # save the dataframe to a new file
        #path_mod = dir_results / Path('230401.1747.air.3.frag.1647s1661s_EI_only_min_4_peaks_per_compound_mod_'+str(factor)+'x_uncertainty.txt')
        path_mod = dir_results / Path('230401.1747.air.3.frag.1592s1600s_'+str(factor)+'x_uncertainty.txt')
        df.to_csv(path_mod, sep='\t', header=True, index=False)
        print(df)
        # replace values in column "mass" with the nearest value in the dictionary

        ## END OF MODIFICATION OF INPUT DATA##


        # run alpinac with all components
        make_identification_from_file(path_file=path_mod, show_plt_figures= False) 
        print("alpinac run wit all possible atoms")
        col_name_results = "alpinac_results"
        suffix = ""


if check_largest_peaks_only == True:
    # get the largest peaks only
    dir_results = Path(r"C:\Users\kaho\Desktop\data\data_Empa\tests")
    #df = df.sort_values(by=['int'], ascending=False)
    #df = df.drop_duplicates(subset=['mass'], keep='first')
    # save the dataframe to a new file
    path_mod = dir_results / Path('230401.1747.air.3.frag.1592s1600s.txt')
    #df.to_csv(path_mod, sep='\t', header=True, index=True)
    print(df)
    # run alpinac with all components
    make_identification_from_file(path_file=path_mod, show_plt_figures= False) 
    print("alpinac run wit all possible atoms")
    col_name_results = "alpinac_results"
    suffix = "_largest_peaks_only"
