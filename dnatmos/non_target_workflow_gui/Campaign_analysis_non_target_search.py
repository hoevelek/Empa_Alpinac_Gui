# Plot a pdf report file 

from email import charset
import io
import logging
from pathlib import Path
import tempfile
from fpdf import FPDF
from PIL import Image
import pandas as pd
import numpy as np
from alpinac_gui.matchms_funcs import formula_from_namestring, formula_to_latex, halocarbname2formula, plot_mass_spectra, load_from_mgf, metadata_processing, peak_processing, AlpinacData, halocarbname2formula
from pyvalem.formula import Formula
from molmass import Formula as chem_formula
import pubchempy as pcp
from matplotlib import pyplot as plt
from matplotlib import axes
from db_reader import AlpinacFragmentsReader, JdxReader, AbstractReader
from matchms.Spectrum import Spectrum
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import chardet

# load the datatable
dir_results = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230401.1747.air.3")
# Read from excel
table_data = pd.read_excel(dir_results / "campaign_analysis_data_vertical_layout_with_colours.xlsx", header=0, index_col=0)
cmp_ids_to_check = [706, 610, 827,851,730,745, 20, 835, 181]


import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

plot_single_spectra = True

# Create a PDF file for the multipanel image
pdf_filename = "multipanel_image.pdf"
pdf = PdfPages(pdf_filename)

if plot_single_spectra:
    ind_plot = 1
    row_ind = range(ind_plot)[0]
    row_indici = [2, 5]

    for row_ind in row_indici:
        # Rest of the code...

        # Create a new figure for each subplot
        fig, axes = plt.subplots(2, 1, figsize=(20, 20))
        ax1, ax2 = axes

        # Plot the first subplot
        plot_mass_spectra(ax1, spectrums[0],
                          db_spectra,
                          no_plt_legend=1)
        
        # Plot the second subplot
        for i in range(len(alp_spec)):
            min_ax = np.nanmin([np.nan if x == 0 else x for x in alp_spec_unidentified[i].peaks.mz])
            ax2.set_xlim([min_ax-1, max(alp_spec[i].peaks.mz)+1])

            alp_spec[i].plot_stem(ax2)
            alp_spec_unidentified[i].plot_stem(ax2, col="orange", label=["meas. unident.", "meas. unident."])
            ax2.set_title("File " + path_alp_compounds_dirs[0].parent.parent.name + ": " + path_alp_compounds_dirs[i].name + " of " + str(len(alp_spec)) + " in RT range " + str(path_alp_spec_RT_ranges[RT_range_ind][0]) + "s to " + str(path_alp_spec_RT_ranges[RT_range_ind][1]) + "s (Alpinac assigned fragments)")

        # Save the figure to the PDF file
        pdf.savefig(fig)
        plt.close(fig)

# Close the PDF file
pdf.close()























plot_single_spectra = True
# SINGLE SPECTRUM PLOTS
if plot_single_spectra:

    ind_plot = 1
    row_ind = range(ind_plot)[0]
    row_indici = [2,5]

    for row_ind in row_indici:

        # define globally
        formulas_ms = ['XX']*3

        # search for libraries
        db_dir = Path(r"G:\503_Themen\Klima\TargetList")
        df_peaks = pd.read_excel(db_dir / "TargetList.xlsx")

        # %% Load the data bases
        dbs: dict[str, dict[str, Spectrum]] = {
        "nist": JdxReader(db_dir / "nist_db").read_all(),
        "myriam": AlpinacFragmentsReader(db_dir / "myriam_ei").read_all(),
        }

        high_score_subs = table_data.iloc[row_ind]['0_matching_cmp'].replace('(n)', 'nist').replace('(e)', 'myriam').split('\n')
        database_specs = [elem.split()[0] for elem in high_score_subs]
        database_orgs = [elem.split()[1] for elem in high_score_subs]
        db_spectra = []
        db_spectra_names = []

        # get the spectra of the 3 highest scoring compounds
        for i in range(3):
            spec = dbs[database_orgs[i]][database_specs[i]]
            db_spectra.append(spec)
            db_spectra_names.append(str(database_specs[i]) + " (" + str(database_orgs[i]) + ")")


        #ind_plot = table_data.shape[0]
        #for row_ind in row_indici:
        pdf.add_page(orientation='L')
        pdf.cell(0, 10, 'Largest compound no '+str(row_ind+1), 0, 1, 'C', 1)
        pdf.ln(10)
        plot_header_and_corresponding_line(row_ind)

        # load mfg file of cmp
        mfg_file = dir_results/Path('matchms_spectra') / Path(str(table_data.iloc[row_ind]['cmp_id'])+'.mgf')
        # spectrums is a Python list of matchms-type Spectrum objects
        spectrums = list(load_from_mgf(str(mfg_file)))
        spectrums[0].metadata['compound_name'] = 'experimental'
        spectrums[0].metadata['formula'] = 'C'

        
        # process meta data, normalize intensities, reject intensities below a limit and select a certain range
        #spectrums = [metadata_processing(s) for s in spectrums]
        spectrums = [peak_processing(s) for s in spectrums]

        db_spectra = [peak_processing(s) for s in db_spectra]

        # get the spectra corresponding to the entries in the datatable
        # plot spectrum
        #spectrums.plot_mass_spectrum(spectrums[0], mirror_intensity=False, ax=axes, peak_color="C0")

        fig, axes = plt.subplots(1, 1)
        ax = axes
        # plot_spectrum(spectrums[0], mirror_intensity=False, ax=axes, peak_color="C0")
        #y_max = ax.get_ylim()[1]

        plot_mass_spectra(ax, spectrums[0], 
                            db_spectra,
                            no_plt_legend = 1)
        plt.show()
        # add logger
        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)

        path_N_yes_no = "*_EI_only_min_4_peaks_per_compound_no_N"
        path_N_yes_no = "*_EI_only_min_4_peaks_per_compound"
        path_alp_spec = dir_results #/ Path(dir_results.stem)
        # list all directories ending with  "_EI_only_min_4_peaks_per_compound_no_N"
        path_alp_spec_dirs = list(path_alp_spec.glob(path_N_yes_no))
        path_alp_spec_RT_ranges = [str(elem).split(".frag.")[1].split("s")[0:2] for elem in path_alp_spec_dirs]
        # convert string ranges to int ranges
        path_alp_spec_RT_ranges = [[int(elem[0]), int(elem[1])] for elem in path_alp_spec_RT_ranges]
        path_alp_spec_RT_mean = [np.mean(elem) for elem in path_alp_spec_RT_ranges]
        # get the RT of the experimental spectrum
        RT = int(table_data.iloc[row_ind]['RT'])
        # get the index of the RT range that is closest to the RT of the experimental spectrum
        RT_range_ind = np.argmin([abs(RT - elem) for elem in path_alp_spec_RT_mean])
        path_alp_compounds = dir_results / Path(dir_results.name + ".frag." + str(path_alp_spec_RT_ranges[RT_range_ind][0]) + "s" + str(path_alp_spec_RT_ranges[RT_range_ind][1]) + "s" +
                                                str(path_N_yes_no).split("*")[1])
        #check if fodler exists
        #if not os.path.exists(path_alp_compounds):
        #    print("Folder does not exist")
        print(path_alp_compounds)
        # list all directories starting with  "Compound":
        path_alp_compounds_dirs = list(path_alp_compounds.glob("Compound_*"))
        # get the index of the RT range that is closest to the RT of the experimental spectrum
        # resets graphics
        plt.close('all')



        if path_alp_compounds_dirs:

            alp_spec = [AlpinacData.from_alpinac(str(elem)+"/results_file_mygu.txt") for elem in path_alp_compounds_dirs]
            alp_spec_unidentified = [AlpinacData.from_alpinac_unidentified(str(elem)+"/results_file_mygu.txt") for elem in path_alp_compounds_dirs]
            alp_input_spec = [AlpinacData.from_alpinac(str(elem.parent)+".txt") for elem in path_alp_compounds_dirs]

            if alp_spec:
                
                fig, ax = plt.subplots(figsize = (20,10))
                for i in range(len(alp_spec)):
                    # all unidentified will be plotted at 0, limited axes to min(alp_spec[i].spectrum.intensities > 0) and max(values)
                    
                    min_ax = np.nanmin([np.nan if x == 0 else x for x in alp_spec_unidentified[i].peaks.mz])
                    ax.set_xlim([min_ax-1, max(alp_spec[i].peaks.mz)+1])

                    alp_spec[i].plot_stem(ax)
                    alp_spec_unidentified[i].plot_stem(ax, col="orange", label = ["meas. unident.", "meas. unident."])
                    # set title of figure
                    ax.set_title("File " + path_alp_compounds_dirs[0].parent.parent.name + ": " + path_alp_compounds_dirs[i].name + " of " + str(len(alp_spec)) + " in RT range " + str(path_alp_spec_RT_ranges[RT_range_ind][0]) + "s to " + str(path_alp_spec_RT_ranges[RT_range_ind][1]) + "s (Alpinac assigned fragments)")

                    # buf1 = io.BytesIO()
                    # canvas = FigureCanvas(fig)
                    # #fig.savefig(buf3, format='png')
                    # canvas.print_figure(buf1, format='png')
                    # buf1.seek(0)
                    # pil1_image = Image.open(buf1)
                add_matplotlib_plot_to_pdf(fig, pdf)

                plt.show()
            
            # # save plot to image object
            # buf = io.BytesIO()
            # canvas = FigureCanvas(fig)
            # fig.savefig(buf, format='png')
            # buf.seek(0)
            #alp_spec[0].plot(ax)
            #alp_spec[0].peaks.mz
            #alp_spec[0].peaks.intensities
            #alp_spec[0].metadata

            # add peaks of corrsponding alpinac input file before running alpinac
            # open csv file using first line as colnames:
            alp_input_raw = pd.read_csv(str(path_alp_compounds_dirs[0].parent)+".txt", sep = "\t", header = 0)
            # get column names
            alp_input_raw.columns
            # get alp_input_raw column which is called 'mass'
            
            # get mz and intensity values
            alp_input_raw['area']
            # do stem plot in grey color

            # get corresponding spectra from matchms_spectra directory
            # get the index of RT column of the excel file which is closest to RT of experimental spectrum
            RT_ind = np.argmin([abs(RT - int(elem)) for elem in table_data['RT']])
            # get cmp_id of excel file which is nearest to RT of experimental spectrum:
            cmp_id = table_data.iloc[RT_ind]['cmp_id']
            # get the path to the mgf file ending with cmp_id from matchms_spectra directory
            path_mgf_files = list(Path(dir_results/"matchms_spectra").glob("*.mgf"))
            # load mgf file ending with cmp_id from matchms_spectra directory
            mgf_file = [elem for elem in path_mgf_files if str(elem).endswith(str(cmp_id)+".mgf")][0]
            # load mgf file as spectrum
            spectrums = list(load_from_mgf(str(mgf_file)))
            # normalize spectrum
            spectrums = [(spectrum) for spectrum in spectrums]
            spectrums = [peak_processing(s) for s in spectrums]



        


            #fig, ax = plt.subplots(1, 1)
            fig, ax = plt.subplots(figsize=(20, 10))
            # set lower xlimit to 1

            plot_mass_spectra(ax, spectrums[0], col="lightgrey", label = ["py peakfinder extr.", "py peakfinder extr."])

            plot_mass_spectra(ax, alp_input_spec[0], 
                            db_spectra,
                            no_plt_legend = 1, label_peaks = "numeric")




            # set title of figure
            db_spectra_str = ""
            for ind in range(len(db_spectra_names)):
                str_add = str(ind) + ": "  + db_spectra_names[ind] + ", "
                db_spectra_str = db_spectra_str + str_add
            # reomve last comma
            db_spectra_str = db_spectra_str[:-2]    

            
            ax.set_title("File " + path_alp_compounds_dirs[0].parent.parent.name + ": " + path_alp_compounds_dirs[i].name + " of " + str(len(alp_spec)) + " in RT range " + str(path_alp_spec_RT_ranges[RT_range_ind][0]) + "s to " + str(path_alp_spec_RT_ranges[RT_range_ind][1]) + "s (Cosine similarity to database spectra: " + db_spectra_str + ")")

            buf2 = io.BytesIO()
            canvas = FigureCanvas(fig)
            #fig.savefig(buf3, format='png')
            canvas.print_figure(buf2, format='png')
            buf2.seek(0)
            pil2_image = Image.open(buf2)
            add_matplotlib_plot_to_pdf(fig, pdf)

            plt.show()



            fig, ax = plt.subplots(figsize=(20, 10))
            plot_mass_spectra(ax, alp_input_spec[0], spectrums[0], no_plt_legend = 1, label_peaks = "numeric", label = ["Alpinac extraction", "py peakfinder extr."])
            #plt.show()

            #plt.figure(figsize=(10, 10))
            #fig, ax = plt.subplots(1, 1, figsize=(20, 10), dpi = 150)
            
            # plot_spectrum(spectrums[0], mirror_intensity=False, ax=axes, peak_color="C0")
            #y_max = ax.get_ylim()[1]
            #alp_spec[0].plot_stem(ax, label_peaks = "numeric")
            ax.set_title("File " + path_alp_compounds_dirs[0].parent.parent.name + ": " + path_alp_compounds_dirs[i].name + " of " + str(len(alp_spec)) + " in RT range " + str(path_alp_spec_RT_ranges[RT_range_ind][0]) + "s to " + str(path_alp_spec_RT_ranges[RT_range_ind][1]) + "s (Alpinac (upper, calibrated) vs peakfinder (lower)")

            add_matplotlib_plot_to_pdf(fig, pdf)

            plt.show()
        
        




        buf3 = io.BytesIO()
        canvas = FigureCanvas(fig)
        #fig.savefig(buf3, format='png')
        canvas.print_figure(buf3, format='png')
        buf3.seek(0)
        pil3_image = Image.open(buf3)


