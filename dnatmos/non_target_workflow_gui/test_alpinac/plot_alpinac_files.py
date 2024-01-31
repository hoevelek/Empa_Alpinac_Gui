from matplotlib import pyplot as plt
from pathlib import Path
import pandas as pd
from alpinac_gui.matchms_funcs import formula_from_namestring, formula_to_latex, halocarbname2formula, plot_mass_spectra, load_from_mgf, metadata_processing, peak_processing, AlpinacData, halocarbname2formula



path_result = r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230401.1747.air.3\230401.1747.air.3.frag.1735s1743s_EI_only_min_4_peaks_per_compound\Compound_2\results_file_mygu.txt"
title = "230401.1747.air.3.frag.1735s1743s, compound 2"

path_result = r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230401.1747.air.3\230401.1747.air.3.frag.1735s1743s_EI_only_min_4_peaks_per_compound\Compound_3\results_file_mygu.txt"
title = "230401.1747.air.3.frag.1735s1743s, compound 3"

alp_spec = [AlpinacData.from_alpinac(path_result)]
fig, ax = plt.subplots()
alp_spec[0].plot_stem(ax)
ax.set_title(title)

fig, ax = plt.subplots()
alp_spec[0].plot_stem(ax, label_peaks="numeric")
ax.set_title(title)