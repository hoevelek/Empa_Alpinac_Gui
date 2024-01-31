import pandas as pd
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
import matplotlib.ticker as ticker


#from dnathmos_programms import *
# TODO import from utils file
def find_block(ind, start_rows, len_table):
    # find block of rows where ind is in
    # start_rows: list of row indices where a new block starts
    # ind: index of row to be found in which block
    # returns: start and end index of block
    for n_start in start_rows:
        try:
            n_end = start_rows[np.where(start_rows > n_start)[0]][0]
        except:
            n_end = len_table
        if ind >= n_start and ind < n_end:
            return n_start, n_end
    return None


# FUNCTIONS

def get_experimental_input_spectra_and_corresponding_databasespectra(table_data_block, dbs, cmp_id):
    high_score_subs = [cs_name.replace('(n)', 'nist').replace('(e)', 'myriam') for cs_name in table_data_block['cs_matching_cmp']]

    #high_score_subs = table_data.iloc[row_ind]['0_matching_cmp'].replace('(n)', 'nist').replace('(e)', 'myriam').split('\n')
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
    #pdf.add_page(orientation='L')
    #pdf.cell(0, 10, 'Largest compound no '+str(row_ind+1), 0, 1, 'C', 1)
    #pdf.ln(10)
    #plot_header_and_corresponding_line(row_ind)

    # load mfg file of cmp
    mfg_file = dir_results/Path('matchms_spectra') / Path(str(cmp_id)+'.mgf')
    # spectrums is a Python list of matchms-type Spectrum objects
    spectrums = list(load_from_mgf(str(mfg_file)))
    spectrums[0].metadata['compound_name'] = 'experimental'
    spectrums[0].metadata['formula'] = 'C'


    # process meta data, normalize intensities, reject intensities below a limit and select a certain range
    #spectrums = [metadata_processing(s) for s in spectrums]
    spectrums = [peak_processing(s) for s in spectrums]

    db_spectra = [peak_processing(s) for s in db_spectra]
    return spectrums[0], db_spectra

# if main

if __name__ == "__main__":

    # PLOT example spectra
    dir_example = r"G:\503_Themen\Klima\ecToFKampagne\Gemessene_Proben\fuer_Michelle_tmp\230310.0724.tank.5_EI_CI_output_678_mod\Compound_0\results_file_mygu.txt"
    dir_example = r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303_second_twin\230310.0724.tank.5\alpinac_results\230310.0724.tank.5_EI_CI_output_101_mod\Compound_3\results_file_mygu.txt"
    dir_example = r"G:\503_Themen\Klima\kaho\results_campaign2023\Givaudon_real_duplicate_2_first_run_discrete_fitting\230310.0724.tank.5\alpinac_results\230310.0724.tank.5_EI_CI_output_631_mod\Compound_3\results_file_mygu.txt"
    file_out = r"C:\Users\kaho\polybox\Paper\TOFWerk\Plots\230310.0724.tank.5_cmp_631.svg"
    # load alpinac data
    alp_spec = [AlpinacData.from_alpinac(dir_example), AlpinacData.from_alpinac(dir_example, ionization_mode = "CI")]

    # get EI spectra
    

    if alp_spec:

        ### the following is just for quick plottin
    
        fig, ax = plt.subplots(figsize = (20,10))
        tick_spacing = 5
        ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

        min_ax = 300
        for i in range(len(alp_spec)):
            # all unidentified will be plotted at 0, limited axes to min(alp_spec[i].spectrum.intensities > 0) and max(values)
            min_ax = np.nanmin([min_ax, (np.nanmin([np.nan if x == 0 else x for x in alp_spec[i].peaks.mz]))])
            #try:
            #    min_ax = np.nanmin(min_ax, np.nanmin([np.nan if x == 0 else x for x in alp_spec_unidentified[i].peaks.mz]))
            #except:
            #    pass
            ax.set_xlim([min_ax-1, max(alp_spec[i].peaks.mz)+1])

            alp_spec[i].plot_stem(ax, col = "C" + str(i))
            #alp_spec[i].plot_stem(ax, col = "C" + str(i), label_peaks = "numeric")

            #try:
                #alp_spec_unidentified[i].plot_stem(ax, col="orange", label = ["meas. unident.", "meas. unident."])
            #except:
            #    pass
            # set title of figure
            #ax.set_title("File " + path_alp_compounds_dirs[0].parent.parent.name + ": " + path_alp_compounds_dirs[i].name + " of " + str(len(alp_spec)) + " in RT range " + str(path_alp_spec_RT_ranges[RT_range_ind][0]) + "s to " + str(path_alp_spec_RT_ranges[RT_range_ind][1]) + "s (Alpinac assigned fragments)")

            # buf1 = io.BytesIO()
            # canvas = FigureCanvas(fig)
            # #fig.savefig(buf3, format='png')
            # canvas.print_figure(buf1, format='png')
            # buf1.seek(0)
            # pil1_image = Image.open(buf1)
        #add_matplotlib_plot_to_pdf(fig, pdf)

        # set tick distance to 1
        

        plt.show()

    # save as inkscape svg
    fig.savefig(file_out, format = "svg")




    # print mass spectrum
    print(alp_spec.peaks.mz)
    print(alp_spec.peaks.intensities)
    # plot alpinac data
    fig, ax = plt.subplots(figsize=(20, 10))
    alp_spec.plot(ax)
    plt.show()





    # PATHS
    dir_results = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230401.1747.air.3")
    file_report = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230401.1747.air.3\230401.1747.air.3.pdf")
    file_excel = Path(dir_results / "230401.1747.air.3_peak_identification_results_mod_kaho.xlsx")
    path_target_list_org = Path(r"G:\503_Themen\Klima\TargetList\TargetList.csv") # for RT
    path_target_list = Path(r"G:\503_Themen\Klima\TargetList\TargetList_extended.csv") # for RT

    # add logger
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    # SINGLE SPECTRUM PLOTS
    #if plot_single_spectra:

    # load table data
    table_data = pd.read_excel(dir_results / "campaign_analysis_data_vertical_layout.xlsx", sheet_name = "Sheet1")

    # load the libraries
    # search for libraries
    db_dir = Path(r"G:\503_Themen\Klima\TargetList")
    #df_peaks = pd.read_excel(db_dir / "TargetList.xlsx")
    df_peaks = pd.read_csv(path_target_list_org)

    # %% Load the data bases
    dbs: dict[str, dict[str, Spectrum]] = {
    "nist": JdxReader(db_dir / "nist_db").read_all(),
    "myriam": AlpinacFragmentsReader(db_dir / "myriam_ei").read_all(),
        }


    ind_plot = 1
    row_ind = range(ind_plot)[0]
    row_indici = [2,5]
    cmp_indici = [706, 610, 827, 851, 730]
    # tab
    # find corresponding row_indici in table_data cmp_id
    row_indici = []
    cmp_ind_dict = {}
    for cmp_ind in cmp_indici:
        row_ind = table_data[table_data['cmp_id'] == cmp_ind].index[0]
        row_indici.append(row_ind)
        cmp_ind_dict[row_ind] = cmp_ind

    for row_ind in row_indici:
        # PLOT CS MATCHING

        start_rows = table_data['cmp_id'][table_data['cmp_id'].isna() == False].index
        # find index where new block starts

        block_ind = find_block(ind = row_ind, start_rows = start_rows, len_table = table_data.shape[0])

        # get block from table_data
        table_data_block = table_data.iloc[block_ind[0]:block_ind[1]]
        # get the spectra corresponding to the entries in the datatable
        # plot spectrum
        #spectrums.plot_mass_spectrum(spectrums[0], mirror_intensity=False, ax=axes, peak_color="C0")
        formulas_ms = table_data_block['formula']

        # get the spectra corresponding to the entries in the datatable
        exp_spectrum, db_spectra = get_experimental_input_spectra_and_corresponding_databasespectra(table_data_block, dbs, cmp_ind_dict[row_ind])


        fig, axes = plt.subplots(1, 1)
        ax = axes
        # plot_spectrum(spectrums[0], mirror_intensity=False, ax=axes, peak_color="C0")
        #y_max = ax.get_ylim()[1]

        plot_mass_spectra(ax, exp_spectrum, 
                        db_spectra,
                        no_plt_legend = 1)
        # add title
        ax.set_title("File " + dir_results.name + ": " + str(cmp_ind_dict[row_ind]) + " (exp. vs " + str(database_specs[0]) + " (" + str(database_orgs[0]) + "), " + str(database_specs[1]) + " (" + str(database_orgs[1]) + "), " + str(database_specs[2]) + " (" + str(database_orgs[2]) + "))")
        plt.show()

        # save plot to image object
        #buf = io.BytesIO()
        #canvas = FigureCanvas(fig)
        #fig.savefig(buf, format='png')
        #buf.seek(0)
        #pil_image = Image.open(buf)



        # PLOT ALPINAC FRAGMENTS

        path_N_yes_no = "*_EI_only_min_4_peaks_per_compound_no_N"
        path_N_yes_no = "*_EI_only_min_4_peaks_per_compound"
        alpinac_mode_folder = "EI"

        #def plot_alpinac_fragments(row_ind, alpinac_mode_folderS):
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
        #plt.close('all')



        if path_alp_compounds_dirs:

            alp_spec = [AlpinacData.from_alpinac(str(elem)+"/results_file_mygu.txt") for elem in path_alp_compounds_dirs]
            alp_spec_unidentified = [AlpinacData.from_alpinac_unidentified(str(elem)+"/results_file_mygu.txt") for elem in path_alp_compounds_dirs]
            alp_input_spec = [AlpinacData.from_alpinac(str(elem.parent)+".txt") for elem in path_alp_compounds_dirs]

            if alp_spec:
                
                fig, ax = plt.subplots(figsize = (20,10))
                for i in range(len(alp_spec)):
                    # all unidentified will be plotted at 0, limited axes to min(alp_spec[i].spectrum.intensities > 0) and max(values)
                    min_ax = np.nanmin([np.nan if x == 0 else x for x in alp_spec[i].peaks.mz])
                    try:
                        min_ax = np.nanmin(min_ax, np.nanmin([np.nan if x == 0 else x for x in alp_spec_unidentified[i].peaks.mz]))
                    except:
                        pass
                    ax.set_xlim([min_ax-1, max(alp_spec[i].peaks.mz)+1])

                    alp_spec[i].plot_stem(ax)
                    try:
                        alp_spec_unidentified[i].plot_stem(ax, col="orange", label = ["meas. unident.", "meas. unident."])
                    except:
                        pass
                    # set title of figure
                    ax.set_title("File " + path_alp_compounds_dirs[0].parent.parent.name + ": " + path_alp_compounds_dirs[i].name + " of " + str(len(alp_spec)) + " in RT range " + str(path_alp_spec_RT_ranges[RT_range_ind][0]) + "s to " + str(path_alp_spec_RT_ranges[RT_range_ind][1]) + "s (Alpinac assigned fragments)")

                    # buf1 = io.BytesIO()
                    # canvas = FigureCanvas(fig)
                    # #fig.savefig(buf3, format='png')
                    # canvas.print_figure(buf1, format='png')
                    # buf1.seek(0)
                    # pil1_image = Image.open(buf1)
                #add_matplotlib_plot_to_pdf(fig, pdf)

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

            # PLOT ALPINAC INPUT    

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
            RT_ind = np.argmin([abs(RT - int(elem)) for elem in table_data['RT'] if np.isnan(elem) == False])
            # get cmp_id of excel file which is nearest to RT of experimental spectrum:
            cmp_id = cmp_ind_dict[row_ind]
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
            #add_matplotlib_plot_to_pdf(fig, pdf)

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

        #add_matplotlib_plot_to_pdf(fig, pdf)

        plt.show()

            




        buf3 = io.BytesIO()
        canvas = FigureCanvas(fig)
        #fig.savefig(buf3, format='png')
        canvas.print_figure(buf3, format='png')
        buf3.seek(0)
        pil3_image = Image.open(buf3)




            # pdf.add_page('P')


            # # add png to pdf

            # with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp_file:
            #     pil3_image.save(tmp_file.name)
            #     tmp_file.seek(0)
            #     image_path = tmp_file.name
                
            #     # Add the image to fpdf
            #     if pdf.get_y() + pil3_image.height > pdf.h:
            #         pdf.add_page('P')

            #     #shrink image if too large
            #     if pil3_image.width > 20000:
            #         pil3_image = pil3_image.resize((200, 100), Image.ANTIALIAS)
                    

            #     # plot image
            #     pdf.image(image_path, x=0, y=0, w=200, h=100)
            #     # close image
            #     tmp_file.close()

            #pdf.image(pil_image, x = 0, y = 50, w = 200, h = 100)
            #pdf.image(pil2_image, x = 0, y = 150, w = 200, h = 100)
            #pdf.add_page()
            #pdf.image(pil3_image, x = 0, y = 250, w = 200, h = 100)
            



            # alp_spec_unidentified = [AlpinacData.from_alpinac_unidentified(str(path_alp_compounds)+"/results_file_mygu.txt") for elem in path_alp_compounds_dirs]
            # if alp_spec_unidentified:
            #     fig, axes = plt.subplots(1, 1)
            #     ax = axes
            #     for i in range(len(alp_spec_unidentified)):
            #         alp_spec_unidentified[i].plot(ax)
            #     plt.show()


    #pdf.output(file_report, 'F')
    #This will generate a report PDF that includes the PNG results. You can modify the code to add more PNG images, text, or other elements to your report.

    #print acutal file path
    #print(file_report)