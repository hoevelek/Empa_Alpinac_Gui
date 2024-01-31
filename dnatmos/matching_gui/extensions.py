# -*- coding: utf-8 -*-
"""
Created on Tue May 17 09:09:01 2022

@author: kaho
"""

from matplotlib import pyplot as plt
import pandas as pd
import matchms as ms
import re
from pathlib import Path
import glob


from pyvalem.formula import Formula

def formula_to_latex(formula:str):
    latex_formula = Formula(formula).latex
    return latex_formula

def formula_from_namestring(namestring: str):
    name_parts = re.split(r'[., \-!?:_]+', namestring)
    name = None
    for piece in name_parts:
        try: 
            name = Formula(piece)
        except:
            continue
    return name

def plot_assigned_alpinac_spectrum(ax, df: pd.DataFrame, substance: str = "unknown"): #TODO df to alpinac class
    mz_pd =df['Exact mass [m/z]']
    mz = mz_pd.to_numpy()
    intensities_pd = df['% of Total Intensity']
    intensities = intensities_pd.to_numpy()/100
    spectrum = ms.Spectrum(mz, intensities)

    '''Plots a first, e.g. the experimental unknown spectra against the a list
    of other spectra to compare (upside down)'''
    markerline, stemline, baseline = ax.stem(
        spectrum.peaks.mz, spectrum.peaks.intensities, use_line_collection = True, label = "theor.")
    plt.setp(stemline, linewidth=1.0,color='C0')
    plt.setp(baseline, linewidth=0.5, color="black")
    plt.setp(markerline, markersize=1.5,color='C0')
    ax.set_xlabel(r"$m\//z$ (Da)")
    ax.set_ylabel(r"Norm. Intensity (arb. u.)")
    ax.set(ylim=(0, 1.2))
    # draw assignments
    for d, l, r in zip(spectrum.peaks.mz, spectrum.peaks.intensities, df['Formula'].to_numpy()):
        try: 
            Formula(r)
        except:
            continue
        ax.annotate(r'${}$'.format(formula_to_latex(r)), xy=(d, l), xytext=(2,0),
                    textcoords="offset points", va = "bottom", ha="center")
    #ax.set_title(r'${}$'.format(formula_to_latex(str(formula_from_namestring(Path(path).name)))))
    ax.set_title(r'${}$'.format(str(substance)))
    # draw experimental lines
    for i, xc in enumerate(df['Meas. mass [m/z]']):
        if i==0: 
            plt.axvline(x=xc, color = 'C1', linewidth=0.5, alpha = 0.7, label= "meas.")
        else:
            plt.axvline(x=xc, color = 'C1', linewidth=0.5, alpha = 0.7)
    ax.legend()

if __name__ == "__main__":
    fig = plt.figure()
    ax = plt.subplot(2,1,1)

    path = r"C:\Users\kaho\polybox\ALPINAC\Alpinac_packing\halosearch\data\nontarget_screening\formulas\200602.1550.std.7.frag.train_C2H6.comp.1275.56_0.txt"

    df = pd.read_csv(path, header = 7,sep = "\t").sort_values(by = 'Exact mass [m/z]')
    f = open(path)
    substance = f.readlines(6)[0].split("_")[-1].split("\n")[0]
    plot_assigned_alpinac_spectrum(ax, df, substance)

    import glob
    files = glob.glob(r"C:\Users\kaho\polybox\ALPINAC\alpinac\data\nontarget_screening\formulas\*")
    print(files)




    fig = plt.figure()
    for i in [1,2]:
        ax = plt.subplot(2,1,i)
        path = files[i]
        df = pd.read_csv(path, header = 7,sep = "\t").sort_values(by = 'Exact mass [m/z]')
        f = open(path)
        substance = f.readlines(6)[0].split("_")[-1].split("\n")[0]
        plot_assigned_alpinac_spectrum(ax, df, substance)

    #df.keys = ['perc_norm_I, exact', 'RT/s', 'Total Intensity', '% of Total Intensity',
    #       'Meas. mass /(m/z)', 'Exact mass/(m/z)', 'Exact mass, u [ppm]',
    #       'Mass diff [ppm]', 'Formula', 'DBE', 'Intensity', 'Intensity, fraction',
    #       'Ionisation mode', 'Spectra name']