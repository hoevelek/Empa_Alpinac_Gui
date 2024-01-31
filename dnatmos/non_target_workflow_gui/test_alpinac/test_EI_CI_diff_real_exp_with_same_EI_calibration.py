import pandas as pd

file = r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230401.1747.air.3\230401.1747.air.3.frag.1682s1690s_EI_CI\Compound_0\results_file_mygu.txt"

#file = r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230401.1747.air.3\230401.1747.air.3.frag.1682s1690s_EI_CI\Compound_0\results_file_mygu.txt"
#file = r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230401.1747.air.3\230401.1747.air.3.frag.2210s2218s_EI_CI\Compound_0\results_file_mygu.txt"
#file = r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230401.1747.air.3\230401.1747.air.3.frag.2210s2218s_EI_CI\Compound_0\results_file_mygu.txt"
# read file as pandas dataframe
columnnames = ["% norm I", "exact	RT [s]",	"Total Intensity",	"% of Total Intensity",	"Meas. mass [m/z]", "Exact mass [m/z]",	"Exact mass u [ppm]", "Mass diff [ppm]", "Formula","DBE", "Intensity",	"Intensity fraction", "Ionisation mode", "Max. allowed adduct",	"Spectra name"]
df = pd.read_csv(file, sep='\t', header=9,  names=columnnames)

# plot column " measured mass" vs. "exact mass" and if Ionisation mode == EI use green color, if Ionisation mode == CI use red color
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1,1)
plt.plot(df["Meas. mass [m/z]"], df["Exact mass [m/z]"], 'x', color='black')
plt.plot(df["Meas. mass [m/z]"][df["Ionisation mode"] == "EI"], df["Exact mass [m/z]"][df["Ionisation mode"] == "EI"], 'x', color='green')
plt.plot(df["Meas. mass [m/z]"][df["Ionisation mode"] == "CI"], df["Exact mass [m/z]"][df["Ionisation mode"] == "CI"], 'x', color='red')
plt.xlabel("measured mass")
plt.ylabel("exact mass")
# add a line with slope 1
plt.plot([0, 200], [0, 200], color='black', linestyle='-', linewidth=1)
plt.show()  

# plot column " Exact mass" vs. "exact mass" and if Ionisation mode == EI use green color, if Ionisation mode == CI use red color
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1,1)
#plt.plot(df["Exact mass [m/z]"], df["Mass diff [ppm]"], 'x', color='black')
plt.plot(df["Exact mass [m/z]"][df["Ionisation mode"] == "EI"], df["Mass diff [ppm]"][df["Ionisation mode"] == "EI"], 'x', color='green')
plt.plot(df["Exact mass [m/z]"][df["Ionisation mode"] == "CI"], df["Mass diff [ppm]"][df["Ionisation mode"] == "CI"], 'x', color='red')
plt.xlabel("Exact mass")
plt.ylabel("Mass diff [ppm]")
# add labels
for i, txt in enumerate(df["Formula"]):
    ax.annotate(txt, (df["Exact mass [m/z]"][i], df["Mass diff [ppm]"][i]))
# add a line with slope 1
plt.plot([0, 200], [0, 0], color='black', linestyle='-', linewidth=1)



# limit df to those fragments which are twice in df[formula]
df = df[df.duplicated(subset=['Formula'], keep=False)]




# plot column " Exact mass" vs. "exact mass" and if Ionisation mode == EI use green color, if Ionisation mode == CI use red color
import matplotlib.pyplot as plt
import numpy as np
fig, ax = plt.subplots(1,1)
df_EI = df[df["Ionisation mode"] == "EI"]
df_CI = df[df["Ionisation mode"] == "CI"]
plt.plot(df_EI["Exact mass [m/z]"], df_EI["Mass diff [ppm]"], 'x', color='green')
plt.plot(df_CI["Exact mass [m/z]"], df_CI["Mass diff [ppm]"], 'x', color='red')
plt.xlabel("Exact mass")
plt.ylabel("Mass diff [ppm]")
# add labels
for i, txt in enumerate(df_EI["Formula"]):
    ax.annotate(txt, (np.array(df_EI["Exact mass [m/z]"])[i], np.array(df_EI["Mass diff [ppm]"])[i]))
for i, txt in enumerate(df_CI["Formula"]):
    ax.annotate(txt, (np.array(df_CI["Exact mass [m/z]"])[i], np.array(df_CI["Mass diff [ppm]"])[i]))

   

# add a line with slope 1
plt.plot([0, 200], [0, 0], color='black', linestyle='-', linewidth=1)



