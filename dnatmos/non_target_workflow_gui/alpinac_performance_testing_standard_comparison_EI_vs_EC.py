from pathlib import Path
import pandas as pd
from molmass import Formula as chem_formula
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# import get_compound_class from alpinac_performance_testing
def get_compound_class(sum_formula, smiles):
    class_compound = "other"
    try:
        # is a CFC if keys "C", and keys "F" or "Cl" are present and other than keys "C","F" or "Cl" there are no other ones
        elements = chem_formula(sum_formula)._elements
    except:
        class_compound = "unknown"
        return class_compound   
    if set(elements.keys()) <= set(["C", "F", "Cl", "H"]): #contains only C, F, Cl, H
        if "C" in elements.keys() : #contains C atoms
            if "H" in elements.keys() : #contains H atoms -> partly halogenated
                if "F" in elements.keys() and "Cl" not in elements.keys():
                    if "=" in smiles:
                        class_compound =  "HFO" # hydrofluoroolefins
                    else:
                        class_compound =  "HFC" # hydrofluorocarbons
                else:
                    class_compound =  "HCFC" # hydrochlorofluorocarbons
            else: # fully halogenated
                # if key "F" is present, and "Cl" is not present, it is a PFC
                if "F" in elements.keys() and "Cl" not in elements.keys(): 
                    class_compound =  "PFC" # perfluorocarbons
                else:
                    class_compound =  "CFC" # chlorofluorocarbons
    if set(elements.keys()) <= set(["C", "H"]): #hydrocarbon
        class_compound = "hydrocarbon"
    # halons
    if set(elements.keys()) <= set(["C", "F", "Cl", "Br", "I", "H"]): #contains only C, F, Cl, H
        if "Cl" in elements.keys() or "Br" in elements.keys() or "I" in elements.keys(): #contains Cl, Br or I atoms
            class_compound =  "halon"

    return class_compound

# load excel

path = Path(r'C:\Users\kaho\polybox\Presentations\Innosuisse')
file = 'reportage.xlsx'
df_EC = pd.read_excel(path / file, sheet_name='EC.txt')
df_EI = pd.read_excel(path / file, sheet_name='EI.txt')
get_compound_class(df_EC['expected_formula'][0], "mm00")

# add column with compound class
df_EC['cmp_class'] = df_EC.apply(lambda x: get_compound_class(x['expected_formula'], "bbb"), axis=1)
df_EI['cmp_class'] = df_EI.apply(lambda x: get_compound_class(x['expected_formula'], "bbb"), axis=1)


# get unique filenames
unique_filenames_EC = df_EC['file'].unique()
unique_filenames_EI = df_EI['file'].unique()

# seperate it into classes
# EC


colors = [(1,1,1), (1,1,1), (1, 0, 0), (0.5, 0, 0), (0, 0.5, 0), (0, 1, 0)]
cmap_user = LinearSegmentedColormap.from_list('RedToGreen', colors, N=6)


from matplotlib.colors import LinearSegmentedColormap

def get_color_map(dict1, dict2, ax, cmp_class=""):
    # Define color mapping
    color_mapping = {
        4: 'darkgreen',
        1: 'darkred',
        2: 'lightcoral',
        3: 'lightgreen'
    }


    # Create a matrix to store colors
    # all keys
    keys_part1 = dict1.keys()
    keys_part2 = dict2.keys()
    keys = list(keys_part1) + list(keys_part2)
    unique_keys = list(set(keys))
    # all values
    matrix = np.zeros((2, len(unique_keys)))



    # Fill the matrix with corresponding colors based on the dictionaries
    for i, key in enumerate(unique_keys):
        #matrix[0, i] = color_mapping.get(dict1[key], 'white')
        #matrix[1, i] = color_mapping.get(dict2[key], 'white')
        if key in keys_part1:
            matrix[0, i] = dict1[key]
        else:
            matrix[0, i] = -1
        if key in keys_part2:
            matrix[1, i] = dict2[key]
        else:
            matrix[1, i] = -1

    # Create the figure and plot
    cax = ax.matshow(matrix, cmap=cmap_user, vmin=-1, vmax=4)

    # Set y-axis labels and positions
    ax.set_yticks([0, 1])
    ax.set_yticklabels(['EI TOF data', 'EC TOF data'])

    # Set x-axis labels vertically
    ax.set_xticks(np.arange(len(unique_keys)))
    ax.set_xticklabels(unique_keys, rotation='vertical')

    # Add a colorbar
    cbar = plt.colorbar(cax)

    # Set labels and title
    plt.xlabel('Standard substances', labelpad=20)
    plt.ylabel('Meas. method')
    #plt.title('Matrix Plot')
    plt.title(cmp_class+"s")
    if cmp_class == "HFC":
        plt.title("HFCs & HFOs")

    # Show the plot
    plt.show()




cmp_classes = df_EC['cmp_class'].unique()





dict_grades_EC_tot = {}
dict_grades_EI_tot = {}

for cmp_class in cmp_classes:
    # get dataframes for each class\
    df_EC_class_cmp = df_EC[df_EC['cmp_class'] == cmp_class]
    df_EI_class_cmp = df_EI[df_EI['cmp_class'] == cmp_class]
    # get unique filenames
    unique_filenames_EC_cmp = df_EC_class_cmp['file'].unique()
    unique_filenames_EI_cmp = df_EI_class_cmp['file'].unique()
    # get the grades for each file and average them
    dict_grades_EC = {}
    dict_grades_EI = {}
    for filename_EI, filename_EC in zip(unique_filenames_EI_cmp, unique_filenames_EC_cmp):
        df_EC_class_cmp_file = df_EC_class_cmp[df_EC_class_cmp['file'] == filename_EC]
        df_EI_class_cmp_file = df_EI_class_cmp[df_EI_class_cmp['file'] == filename_EI]
        # get formulas
        formula_EC = df_EC_class_cmp_file['expected_formula'].iloc[0]
        formula_EI = df_EI_class_cmp_file['expected_formula'].iloc[0]
        # get intersection of formulas
        formula_intersection = set(formula_EC) & set(formula_EI)
        # for each element in intersection get the grades and average them
        grades_EC = []
        grades_EI = []
        for element in formula_intersection:
            # get grades
            grade_EC = df_EC_class_cmp_file[df_EC_class_cmp_file['expected_formula'].str.contains(element)]['grade'].iloc[0]
            grade_EI = df_EI_class_cmp_file[df_EI_class_cmp_file['expected_formula'].str.contains(element)]['grade'].iloc[0]
            # append to list
            grades_EC.append(grade_EC)
            grades_EI.append(grade_EI)
        # average the grades
        mean_grade_EC = sum(grades_EC) / len(grades_EC)
        mean_grade_EI = sum(grades_EI) / len(grades_EI)
        # append to dict
        dict_grades_EC[filename_EC.replace("_EC.txt","").replace("frag_", "")] = mean_grade_EC
        dict_grades_EI[filename_EI.replace("_EI.txt","").replace("frag_", "")] = mean_grade_EI
    # plot file
    dict_grades_EC_tot[cmp_class] = dict_grades_EC
    dict_grades_EI_tot[cmp_class] = dict_grades_EI



    fig, ax = plt.subplots(figsize=(8, 6))
    get_color_map(dict_grades_EI, dict_grades_EC, ax, cmp_class)
    # save figure
    path_save = Path(r'C:\Users\kaho\polybox\Presentations\Innosuisse')
    fig.savefig("matrix_plot_"+cmp_class+".png")
    # set title




############################################################################################################
# PIE PLOTS
############################################################################################################

# flatten dict
dict_grades_EC_tot_flat = {}
dict_grades_EI_tot_flat = {}
for cmp_class in cmp_classes:
    dict_grades_EC_tot_flat.update(dict_grades_EC_tot[cmp_class])
    dict_grades_EI_tot_flat.update(dict_grades_EI_tot[cmp_class])

fig, ax2 = plt.subplots(figsize=(8, 6))
get_color_map(dict_grades_EI_tot_flat, dict_grades_EC_tot_flat, ax2, "All classes")
# save figure
path_save = Path(r'C:\Users\kaho\polybox\Presentations\Innosuisse')
fig.savefig("matrix_plot_all_classes.png")


# calculate how many samples are in each class
# EC
dict_samples_EC = {}
for cmp_class in cmp_classes:
    dict_samples_EC[cmp_class] = len(dict_grades_EC_tot[cmp_class])
# EI
dict_samples_EI = {}
for cmp_class in cmp_classes:
    dict_samples_EI[cmp_class] = len(dict_grades_EI_tot[cmp_class])

# calculate how many samples are in each grade
# EC
dict_grades_EC = {}
for cmp_class in cmp_classes:
    dict_grades_EC[cmp_class] = {}
    for grade in range(1,5):
        dict_grades_EC[cmp_class][grade] = 0
    for grade in dict_grades_EC_tot[cmp_class].values():
        dict_grades_EC[cmp_class][grade] += 1
# EI
dict_grades_EI = {}
for cmp_class in cmp_classes:
    dict_grades_EI[cmp_class] = {}
    for grade in range(1,5):
        dict_grades_EI[cmp_class][grade] = 0
    for grade in dict_grades_EI_tot[cmp_class].values():
        dict_grades_EI[cmp_class][grade] += 1   



# forevery class get the percentages of each grade
dict_grades_EC_perc_tot = {}
dict_grades_EI_perc_tot = {}
# initilaize dict by setting all values to 0 for grades 1,2,3,4
dict_grades_EC_perc_tot[1] = 0
dict_grades_EC_perc_tot[2] = 0
dict_grades_EC_perc_tot[3] = 0
dict_grades_EC_perc_tot[4] = 0
dict_grades_EI_perc_tot[1] = 0
dict_grades_EI_perc_tot[2] = 0
dict_grades_EI_perc_tot[3] = 0
dict_grades_EI_perc_tot[4] = 0

for cmp_class in cmp_classes:
    for grade in range(1,5):
        # add to dict and sum if key already exists
        if grade in dict_grades_EC_perc_tot.keys():
            dict_grades_EC_perc_tot[grade] += dict_grades_EC[cmp_class][grade]
            dict_grades_EI_perc_tot[grade] += dict_grades_EI[cmp_class][grade]

    

# get percentages
# EC
dict_grades_EC_perc = {}
for cmp_class in cmp_classes:
    dict_grades_EC_perc[cmp_class] = {}
    for grade in range(1,5):
        dict_grades_EC_perc[cmp_class][grade] = dict_grades_EC[cmp_class][grade] / dict_samples_EC[cmp_class]   
# EI
dict_grades_EI_perc = {}
for cmp_class in cmp_classes:
    dict_grades_EI_perc[cmp_class] = {}
    for grade in range(1,5):
        dict_grades_EI_perc[cmp_class][grade] = dict_grades_EI[cmp_class][grade] / dict_samples_EI[cmp_class]

# do pie plot for each class with the percentages of each grade in 4x2 subplots
# EC
fig, axs = plt.subplots(2, 5, figsize=(15, 10))
for i, cmp_class in enumerate(cmp_classes):
    ax = axs[1//2, i%5]
    ax.pie(dict_grades_EC_perc[cmp_class].values(), labels=dict_grades_EC_perc[cmp_class].values(), colors=cmap_user(np.arange(6))[2:])
    ax.set_title(cmp_class + ", EC, (" + str(len(dict_grades_EC_tot[cmp_class])) + " sample(s))")   

# plot

ax = axs[1//2, 4]
labels = dict_grades_EC_perc_tot.values()
rounded_labels = [round(elem/sum(labels) , 2)for elem in labels]
ax.pie(dict_grades_EC_perc_tot.values(),labels = rounded_labels, colors=cmap_user(np.arange(6))[2:])
ax.set_title("All classes, EC, (" + str(len(dict_grades_EI_tot_flat)) + " sample(s))")
for i, cmp_class in enumerate(cmp_classes):
    ax = axs[2//2, i%5]
    labels = dict_grades_EI_perc[cmp_class].values()
    rounded_labels = [round(elem , 2) for elem in labels]
    ax.pie(dict_grades_EI_perc[cmp_class].values(), labels=rounded_labels, colors=cmap_user(np.arange(6))[2:])
    ax.set_title(cmp_class + ", EI, (" + str(len(dict_grades_EI_tot[cmp_class])) + " sample(s))")   
ax = axs[2//2, 4]
labels = dict_grades_EI_perc_tot.values()
rounded_labels = [round(elem/sum(labels), 2) for elem in labels]
ax.pie(dict_grades_EI_perc_tot.values(),labels = rounded_labels, colors=cmap_user(np.arange(6))[2:])
ax.set_title("All classes, EI, (" + str(len(dict_grades_EI_tot_flat)) + " sample(s))")

# save figure
path_save = Path(r'C:\Users\kaho\polybox\Presentations\Innosuisse')
fig.savefig("pie_plots_EC_vs_EI.png")

# # get total percentages
# # EC
# dict_grades_EC_perc_tot = {}
# for grade in range(1,5):
#     dict_grades_EC_perc_tot[grade] = 0
# for cmp_class in cmp_classes:
#     for grade in range(1,5):
#         dict_grades_EC_perc_tot[grade] += dict_grades_EC_perc[cmp_class][grade]
# # EI
# dict_grades_EI_perc_tot = {}
# for grade in range(1,5):
#     dict_grades_EI_perc_tot[grade] = 0
# for cmp_class in cmp_classes:
#     for grade in range(1,5):
#         dict_grades_EI_perc_tot[grade] += dict_grades_EI_perc[cmp_class][grade]

# # plot pie chart
# # EC
# fig, ax = plt.subplots()
# ax.pie(dict_grades_EC_perc_tot.values(), labels=dict_grades_EC_perc_tot.values(), colors=['darkred', 'lightcoral', 'lightgreen', 'darkgreen'])
# ax.set_title("EC, (" + str(len(dict_grades_EC_tot_flat)) + " sample(s))")
# # EI
# fig, ax = plt.subplots()
# ax.pie(dict_grades_EI_perc_tot.values(), labels=dict_grades_EI_perc_tot.keys())
# ax.set_title("EI, (" + str(len(dict_grades_EI_tot_flat)) + " sample(s))")

# # set title
# plt.title("All classes")

    # do a matrix plot, x axis are the keys of the dict, y axis are the values of the dict (y=0 is EC, y=1 is EI). If grade is 4 plot it dark green, if grade is 1 plot it dark red, if grade is 2 plot it light red, if grade is 3 plot it light green
    # plot EC
    #fig, ax = plt.subplots()
    #ax.scatter(dict_grades_EC.keys(), dict_grades_EC.values())
    # title is cmp_class
    #ax.set_title(cmp_class + ", EC, (" + str(len(dict_grades_EC)) + " sample(s))")






    
    # add subtitle with sum_formula

        # plot per
        # get grade using mean and std for each

        # append to list


    # do pie plot for each class with the percentages of each grade
    #fig, ax = plt.subplots()
    #ax.pie(dict_grades_EC.values(), labels=dict_grades_EC.values())
    # title is cmp_class
    #ax.set_title(cmp_class + ", EC, (" + str(len(dict_grades_EC)) + " sample(s))")
    # add subtitle with sum_formula


    # plot per
        



    # get grade using mean and std for each 