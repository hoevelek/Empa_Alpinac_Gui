# import cirpy
from pathlib import Path
import cirpy
import pandas as pd
import pubchempy as pcp
import numpy as np
import matplotlib.pyplot as plt
#  

# get a resultsfile
# get a file with the measured data
path_to_file = Path(r"g:\503_Themen\Klima\kaho\results_campaign2023\BASF_single_1\230313.1808.tank.3\230313.1808.tank.3_peak_identification_results_extended.xlsx")
# load file
df_res = pd.read_excel(path_to_file)

# read in the standard list
#std_cmps_names = ["HFC-41", "HFC-1234yf", "HFC-1233zdE", "HFC-1234zeE", PFC-216, c-C4F8O, H-1202, HFC-161, H-2311. And NOVEC-4710, HCFC-151a, HFO-1225ye(E), HFO-1225ye(Z), HFO-1335mzz(Z), HFO-1336mzz(E), HBFO-1233xfB, C5-ketone (potential HFC-227ea contamination), i-C4F10, i-C6F14 (CAS 354-96-1), iC6F14 (CAS 355-04-4).
std_cmps_names = ["HFC-41", "HFC-1234yf", "HFC-1233zdE", "HFC-1234zeE", "PFC-216", "c-C4F8O", "H-1202", "HFC-161", "H-2311", "NOVEC-4710", "HCFC-151a", "HFO-1225ye(E)", "HFO-1225ye(Z)", "HFO-1335mzz(Z)", "HFO-1336mzz(E)", "HBFO-1233xfB", "3-pentanone", "2-pentanone", "i-C4F10", "i-C6F14", "iC6F14"]
# is it really HFC_1234yf? not HFO-1234yf?

# use cirpy resolve to get the cas numbers and inchi keys
# load additional data
file_dir = Path(r"C:\Users\kaho\Desktop\data\data_Empa\databases_Empa\TargetList\TargetList_Additional_information_mod_kaho.xlsx")
df_add = pd.read_excel(file_dir)
# search for std_cmps_names in df_add and get the cas numbers and inchi keys
cas_nos = {}
#    cirpy.resolve('1,1,1,4,4,4-hexafluoro-2-butene', 'cas')

cas_nos['HFC-1234yf']  = '754-12-1' # for HFO (same as HFC)
cas_nos['HFC-1233zdE'] = '102687-65-0'
cas_nos['HFC-1234zeE'] = '29118-24-9'
cas_nos['c-C4F8O'] = '773-14-8'
cas_nos['HFO-1225ye(E)'] = '5595-10-8'
cas_nos['HFO-1225ye(Z)'] = '2252-83-7'
#HFO-1335mzz(Z) #? can not find 1335mzz(Z) in the list is it 1336mzz(Z)?
cas_nos['HFO-1336mzz(E)'] = '66711-86-2'
cas_nos['HBFO-1233xfB'] = '1514-82-5'
#C5-ketone: are this all?
cas_nos['2-pentanone'] = '107-87-9'
cas_nos['3-pentanone'] = '96-22-0'
cas_nos['i-C6F14'] = '354-96-1'
cas_nos['iC6F14'] = '355-04-4'

for cmp in std_cmps_names:
    if cmp not in cas_nos.keys():
        try:
            cas_no = df_add[df_add["Compound"] == cmp]["CAS"].values[0]
        #inchikey = df_add[df_add["Compound"] == cmp]["InChIKey"].values[0]
            cas_nos[cmp] = cas_no
        except IndexError:
            print(cmp)
            cas_nos[cmp] = "unknown"


# try to find cas in df_add and get the RT
RTs = {}
for cmp in cas_nos.keys():
    try:
        RT = df_add[df_add["CAS"] == cas_nos[cmp]]["RT"].values[0]
        RTs[cmp] = RT
    except IndexError:
        print(cmp)
        RTs[cmp] = "unknown"



# get all cmps from df_res
cmps = df_add["Compound"].values

RTs = {}
for cmp in cmps:
    try:
        RT = df_add[df_add["Compound"] == cmp]["RT"].values[0]
        RTs[cmp] = RT
    except IndexError:
        print(cmp)
        RTs[cmp] = "unknown"

# for each class print boiling points vs RTs
# get the boiling points from df_add
bpts = {}
for cmp in cmps:
    try:
        bpt = df_add[df_add["Compound"] == cmp]["Boiling point"].values[0]
        bpts[cmp] = bpt
    except IndexError:
        print(cmp)
        bpts[cmp] = np.nan

# plot the RTs vs bpts for each class
# get the classes
classes = {}
for cmp in cmps :
    try:
        cls = df_add[df_add["Compound"] == cmp]["cmp_class"].values[0]
        classes[cmp] = cls
    except IndexError:
        print(cmp)
        classes[cmp] = np.nan
        # 
# define colours for each class using first python colors
colour_class = {}   
for each_class in set(classes.values()):
    colour_class[each_class] = "C" + str(list(set(classes.values())).index(each_class))
# plot the RTs vs bpts for each class 
fig, ax = plt.subplots(1,1)
# plot the RTs vs bpts for each class
for each_class in set(classes.values()):
    # get the cmps of the class
    cmps = []
    for cmp in classes.keys():
        if classes[cmp] == each_class:
            cmps.append(cmp)
    # get the RTs
    RTs_class = {}
    for cmp in cmps:
        RTs_class[cmp] = RTs[cmp]
    # get the bpts
    bpts_class = {}
    for cmp in cmps:
        bpts_class[cmp] = bpts[cmp]
    # plot
    # remove '?' from bpts and RTs
    # convert to np.array
    bpts_class = np.array(list(bpts_class.values()))
    RTs_class = np.array(list(RTs_class.values()))
    # check if is number


    is_numeric = np.char.isnumeric(bpts_class.astype(str))
    # Check if elements contain a single decimal point
    has_decimal_point = np.char.count(bpts_class.astype(str), '.') == 1
    # Combine the two checks to find numbers with or without a decimal point
    has_minus_sign = np.char.count(bpts_class.astype(str), '-') == 1
    mask1 = is_numeric | has_decimal_point | has_minus_sign
    # check if is number
    is_numeric = np.char.isnumeric(RTs_class.astype(str))
    # check if_number and larger than 0
    is_positiv = is_numeric & (RTs_class.astype(float) > 10)
    # Check if elements contain a single decimal point
    has_decimal_point = np.char.count(RTs_class.astype(str), '.') == 1
    # Combine the two checks to find numbers with or without a decimal point
    has_minus_sign = np.char.count(RTs_class.astype(str), '-') == 1
    mask2 = (is_numeric | has_decimal_point) & ~has_minus_sign
    # apply mask

    bpts_class = bpts_class[mask1 & mask2]
    RTs_class = RTs_class[mask1 & mask2]

    # scatter plot on same y-axis

    ax.scatter(RTs_class.astype(float), bpts_class.astype(float), label=each_class, color = colour_class[each_class])
    # add fit for each class but "other"
    if classes[cmp] != "other":
        # fit a line
        z = np.polyfit(RTs_class.astype(float), bpts_class.astype(float), 3)
        p = np.poly1d(z)
        ax.plot(RTs_class.astype(float),p(RTs_class.astype(float)), color = colour_class[each_class], label = each_class + " fit")
# add legend
plt.legend()

plt.title("TR-bpt behaviour by class")
plt.xlabel("RT (s)")
plt.ylabel("bpt (deg C)")
plt.show()






    #inchikeys.append(inchikey)







std_cmps_cas = ['593-53-3']
std_cmps_inchi = ["InChI=1S/CH3F/c1-2/h1H3"]
for cmp in std_cmps_names[1:]:
    std_cas = cirpy.resolve(cmp, 'cas')
    print(cmp, std_cas)
    std_cmps_cas.append(std_cas)
    inchikey = cirpy.resolve(cmp, "stdinchikey")
    print(cmp, inchikey)
    std_cmps_inchi.append(inchikey)

cas_no = {}
# name to compund
for identifier in std_cmps_names:
    cmps = pcp.get_compounds(identifier, 'name')
    print(identifier, cmps)
    inchis = [cmp.inchi for cmp in cmps]
    cas_no[identifier] = [cirpy.resolve(inchi, 'cas') for inchi in inchis] 


#path_to_std_docx = 