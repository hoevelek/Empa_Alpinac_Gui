




# excecute if called as main
from datetime import datetime
from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from data_analysis_utils import formula_by_name, get_calibrant_peaks, get_mass_spec_for_RT_range, get_ranges_for_summing, get_sum, index_to_mass, load_data, get_RT_spec_for_mass_range
from alpinac.utils_data_extraction import f_tofidx_to_mass, f_mass_to_tofidx
from scipy.interpolate import griddata
from scipy.interpolate import interp1d


def get_GCWorks_target_lines(file_path, meas_datetimes):
    # load data
    df = pd.read_csv(
        file_path,
        # Separator
        sep="\s+",
        # use the line number 1 and 2 as headers # Line 0 is skipped
        header=[1, 2],
        # Ensure the date an time are read as str
        dtype={("-", "date"): str, ("-", "time"): str},
    )
    df[('-', 'datetime')] = (
        pd.to_datetime(
            # Merge the str for date and time
            df["-", "date"] + df["-", "time"],
            format="%y%m%d%H%M%S",
        )
    )

    # get all substances that are in the targetlist
    substances = set([cmp for cmp, var in df.columns if cmp != '-'])
    # get index by checking which measurement datetime corresponds to the datetime in the GCWorks file
    index_meas = [i for i, date in enumerate(df[('-', 'datetime')]) if date in meas_datetimes]
    returnlist = []
    for index in index_meas:
        serie = df.iloc[index]
        # plot for each substance the RT 
        out = {}
        for col in df.columns:
            if col[0] == '-':
                print(col)
                continue

            rt = serie[col]
            out[rt] = col

            df_peaklist = pd.DataFrame(out).T.sort_index()
            # add column to retunr list
            returnlist.append(df_peaklist)

    return returnlist

if __name__ == "__main__":
    # all integer masses from 1 to 300
    masses = np.arange(1, 300)
    #get_chemical_possible_combinations_for_elem_list(['C', 'Cl', 'F'], masses, radicals_allowed=False, path_to_files=r"C:\Users\kaho\polybox\Measurement_Campaign\formulas_non_target", recalculate=True)
    #get_chemical_possible_combinations_for_elem_list(['C', 'Cl', 'F', 'Br', 'H'], masses, radicals_allowed=False, path_to_files=r"C:\Users\kaho\polybox\Measurement_Campaign\formulas_non_target", recalculate=True)



    formula_by_name('463-49-0',{})
    # LOAD DATA
    file_path = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230401.1747.air.3.h5")
    file_path = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303_interesting_file\230311.0125.tank.1.h5")
    file_path = Path(r"G:\503_Themen\Klima\kaho\results_campaign2023\Givaudon_real_duplicate_2_first_run_discrete_fitting\230310.0724.tank.5.h5")

    ei_seg = (3,4)
    ci_seg = (1,2)
    ionisation_dict = {"EI": ei_seg,"CI": ci_seg}
    # These are some masses that we don't want as they have maybe a column bleed
    masses_to_ignore = [185, 187]

    # the idea is to sum up the mass cailbrated data of several files (if it is contained in the data) and extract peaks afterwards
    # find all hdf5 files in folder
    folder = Path(r"C:\Users\kaho\Desktop\data\data_Empa\results_campaign2023")
    all_files = [str(file) for file in folder.glob('**/*.h5')]
    all_files_where_substance_was_found = all_files#
    RT_center_where_substance_was_found = [1898, 1839, 1852, 1995, 2005, 1978, 1881, 1898, 1903, 1538]
    RT_rel_range_where_substance_was_found = [-20, 25]
    folder_to_save_intermediate_results = r"C:\Users\kaho\Desktop\data\data_Empa\unknown_substance\intermediate_results"
    s_per_bin = float(1/5) # 5Hz

    # create folder to save intermediate results if it does not exist
    Path(folder_to_save_intermediate_results).mkdir(parents=True, exist_ok=True)

    # load calibrant data
    calibrant_data_path = r"c:\Users\kaho\Desktop\data\data_Empa\combinatorial_data\formulas_CF_radicals_True.txt"
    calibrant_data = pd.read_csv(calibrant_data_path, sep="\t", header=None, names=["formula", "mass"])
    lines = calibrant_data['mass']
    label = calibrant_data['formula']

    calculate_total_ion_current = False
    if calculate_total_ion_current:
        # calcualate RT spectra
        for i, file in enumerate(all_files_where_substance_was_found):

            for spectrum_type in ['EI', 'CI']:
                data_ionimode, reader = load_data(file, ionisation_dict, spectrum_type)
                # get mass axis
                mass_params = reader.mass_calibration_data[spectrum_type].get_parameters_for_rt(1800)
                mass_axis = f_tofidx_to_mass(np.arange(data_ionimode.shape[1]), mass_params)
                # get sum over RT +- RT_rel_range_where_substance_was_found
                sum_over_mass = get_sum(data_ionimode, sum_type = 'RT', tol=0.1)
                ranges_start = [20,34,36,45,46,47,59,60]
                # convert to mass 
                ranges_start_tofidx = [f_mass_to_tofidx(mass, mass_params) for mass in ranges_start]
                sum_over_mass_halocarbs = get_sum(data_ionimode, sum_type = 'RT', tol=5, range_starts=np.array(ranges_start_tofidx)-5, range_ends=np.array(ranges_start_tofidx)+5)
                np.shape(sum_over_mass_halocarbs)
                RT_axis = np.arange(data_ionimode.shape[0])*0.2
                np.shape(RT_axis)
                # reshape RT_axis to (x, 1)
                RT_axis = RT_axis.reshape((RT_axis.shape[0], 1))
                # reshape sum_over_mass to (x, 1)
                sum_over_mass = sum_over_mass.reshape((sum_over_mass.shape[0], 1))

                # save RT_axis, and sum over mass halocarbons
                data_to_save_RT = np.hstack((RT_axis, sum_over_mass, sum_over_mass_halocarbs))
                # save as txt file
                np.savetxt(folder_to_save_intermediate_results + r"\RT_axis_vs_mass_sum_halocarbs_{}_{}.txt".format(Path(file).stem, spectrum_type), data_to_save_RT, delimiter="\t")

    # load targetlist
    targetlist_path = r"c:\Users\kaho\Desktop\data\data_Empa\databases_Empa\TargetList\TargetList_extended.csv"
    targetlist = pd.read_csv(targetlist_path)
    lines = targetlist['RT']
    labels = targetlist['Substance']

    # load GC work target list
        # load GCWorks target lines
    import pandas as pd
    file_path = r"C:\Users\kaho\Desktop\data\data_Empa\databases_Empa\TargetList\peaklist_ectof.dat"

    # get first two lines
    with open(file_path, 'r') as file:
        # read first two lines
        line1 = file.readline()
        line2 = file.readline()
        line3 = file.readline()
        # get column names
        column_names = line2.split("\t")
        column_names_part2 = line3.split("\t")
        # if there is a "_" in column names part 2 add it at correponding position in column names
        # remove newline character
        #column_names[-1] = column_names[-1].replace("\n", "")
        # get number of columns
        number_of_columns = len(column_names)
        # if there is a "Mr " or "Mrs " in the column names


        # get number of lines
        number_of_lines = sum(1 for line in open(file_path, 'r'))
        # get number of measurements
        number_of_measurements = int(number_of_lines/number_of_columns)
        # create empty dataframe
        GCWorks_target_lines = pd.DataFrame(index = range(number_of_measurements), columns = column_names)
        # read all lines and store in dataframe
        for i, line in enumerate(file.readlines()):
            # remove newline character
            line = line.replace("\n", "")
            # split line by tab
            line_split = line.split("\t")
            # store in dataframe
            GCWorks_target_lines.iloc[int(i/number_of_columns), i%number_of_columns] = line_split[0]


    
    # readlines and replace all "Mr " with "Mr_" and all "Mrs " with "Mrs_" to avoid problems with pandas
    def replace_mr_and_mrs(line):
        line = line.replace("Mr ", "Mr_")
        line = line.replace("Mrs ", "Mrs_")
        line = line.replace("Ms", "Ms")
        line = line.replace("Mr\t", "Mr_")
        line = line.replace("Mrs\t", "Mrs_")
        line = line.replace("Ms\t", "Ms_")
        #
    dataframe = [replace_mr_and_mrs(line) for line in open(file_path, 'r').readlines()]
    with open(file_path, 'r') as file:
        # read every line and replace arbitrary number of spaces with ","
        line = file.read()
        print(line)
        filedata = line.replace("Mr ", "Mr_")
        filedata = line.replace("Mrs ", "Mrs_")
        filedata = line.replace("Ms", "Ms")
        #

    


        



    # Write the file out again
    file_path_formatted = file_path.replace(".dat", "_formatted.dat")
    with open(file_path_formatted, 'w') as file:
        file.write(filedata)
    # read formatted file with pandas
    GCWorks_target_lines = pd.read_csv(file_path_formatted)




    meas_datetime = [Path(names).name for names in all_files_where_substance_was_found]
    # take only first two "." separated parts
    meas_datetime = [date.split(".")[0] + "." + date.split(".")[1] for date in meas_datetime]
    # convert to datetime object
    meas_datetimes = [datetime.strptime(date, '%y%m%d.%H%M') for date in meas_datetime] 





    dfp = get_GCWorks_target_lines(file_path, [meas_datetimes[0]])




    # load files 
    all_files = [str(file) for file in Path(folder_to_save_intermediate_results).glob('**/*.txt')]
    all_files_RT_axis_vs_mass_sum_halocarbs = [file for file in all_files if "RT_axis_vs_mass_sum_halocarbs" in file]
    # only those with EI
    all_files_RT_axis_vs_mass_sum_halocarbs = [file for file in all_files_RT_axis_vs_mass_sum_halocarbs if "EI" in file]

    # gererate list of first 10 python colors
    import matplotlib
    colors = matplotlib.colors.TABLEAU_COLORS
    colors = list(colors.keys())


    import numpy as np
    import scipy.optimize as opt

    
    



    # calibrate RT axis
    calibrating_dict = {}
    calibrating_dict['Benzene'] = 1908.479381
    calibrating_dict['Toluene'] = 2215.94723
    calibrating_dict["Trichlorofluoromethane"]= 1597.372393
    calibrating_dict["Dichlorodifluoromethane"] = 1468.148516
    calibrating_dict["Tetrachloromethane"] = 1740.769826





    exp_calibrant_for_file = {}
    TIC_dict = {}
    for i, file in enumerate(all_files_RT_axis_vs_mass_sum_halocarbs):
        # x is mass axis of first file
        data_txt = np.loadtxt(file, delimiter="\t")
        ncols = np.shape(data_txt)[1]
        RT_axis_eval = data_txt[:,0].reshape((data_txt.shape[0], 1)) #x
        TIC = data_txt[:,1].reshape((data_txt.shape[0], 1)) #y
        # experimental calibrant RT
        exp_calibrant_RT = {}
        # find next max position in range +- 2 in TIC around each calibrant
        for calibrant in calibrating_dict.keys():
            calibrant_RT = calibrating_dict[calibrant]
            # find index of calibrant
            calibrant_index = np.argmin(np.abs(RT_axis_eval - calibrant_RT))
            # find max in TIC +- 2
            max_index = np.argmax(TIC[calibrant_index-10:calibrant_index+10])
            # get RT of max
            max_RT = RT_axis_eval[calibrant_index-10:calibrant_index+10][max_index]
            # calibrate RT axis
            #RT_axis_eval = RT_axis_eval - max_RT + calibrant_RT
            print(max_RT)
            # add to experimental calibrant RT
            exp_calibrant_RT[calibrant] = max_RT[0]
        # add to dict
        exp_calibrant_for_file[Path(file).name] = exp_calibrant_RT
    calibration_param_dict = {}
    # find a correction function a*x + b to match each file with the calibranting dict values
    fig, ax = plt.subplots(1,1)

    for file in exp_calibrant_for_file.keys():
        # get experimental calibrant RT
        exp_calibrant_RT = exp_calibrant_for_file[Path(file).name]
        # get calibrant RT
        calibrant_RT = calibrating_dict
        # do optimization, such that sum of squares of differences from exp. values to theoretical calibrant is minimized
        optim_func = lambda x: sum([(calibrant_RT[calibrant] - (x[0]*exp_calibrant_RT[calibrant] + x[1]))**2 for calibrant in exp_calibrant_RT.keys()])
        from scipy.optimize import minimize
        res = minimize(optim_func, [1, 0])
        a = res.x[0]
        b = res.x[1]
        # store values to calibration_param_dict
        calibration_param_dict[Path(file).name] = [a, b]
        # get calibrated RT axis
        print(exp_calibrant_RT.values())
        print([elem*a + b for elem in exp_calibrant_RT.values()])
        # # do linear fit to get a and b
        # a, b = np.polyfit(list(exp_calibrant_RT.values()), list(calibrant_RT.values()), 1)
        # # store values to calibration_param_dict
        
        # calibration_param_dict[Path(file).name] = [a, b]    
        # # get calibrated RT axis
        # print(exp_calibrant_RT.values())
        # print([elem*a + b for elem in exp_calibrant_RT.values()])
        # # plot against calibrant RT
        # ax.plot(exp_calibrant_RT.values(), [elem*a + b for elem in exp_calibrant_RT.values()], "o")



        exp_calibrant_RT = exp_calibrant_for_file[Path(file).name]
        # get calibrant RT
        calibrant_RT = calibrating_dict
        # do optimization, such that sum of squares of differences from exp. values to theoretical calibrant is minimized
        optim_func = lambda x: sum([(calibrant_RT[calibrant] - (x[0]*exp_calibrant_RT[calibrant] + x[1]))**2 for calibrant in exp_calibrant_RT.keys()])
        from scipy.optimize import minimize
        res = minimize(optim_func, [1, 0])
        a = res.x[0]
        b = res.x[1]













    fig, ax = plt.subplots(1,1)
    for i, file in enumerate(all_files_RT_axis_vs_mass_sum_halocarbs):
        # x is mass axis of first file
        data_txt = np.loadtxt(file, delimiter="\t")
        ncols = np.shape(data_txt)[1]
        RT_axis_eval = data_txt[:,0].reshape((data_txt.shape[0], 1)) #x
        # calibrate RT axis
        a, b = calibration_param_dict[Path(file).name]
        RT_axis_eval = RT_axis_eval*a + b
        TIC = data_txt[:,1].reshape((data_txt.shape[0], 1)) #y
        # other columns are sum over mass halocarbs arrays
        sum_over_mass_halocarbs = data_txt[:,2:ncols]
        # plot
        ax.plot(RT_axis_eval, TIC, color = colors[i])
        ax.plot(RT_axis_eval, -sum_over_mass_halocarbs.sum(axis = 1)*200, linestyle = "dashed", color = colors[i])
        df_peaklist = get_GCWorks_target_lines(file_path, meas_datetimes[i])
        # add lines and labels
        for i, line in enumerate(df_peaklist.index):
            ax.axvline(line, color=colors[i], linestyle="--", linewidth=0.5)
            ax.text(line, 0, df_peaklist.iloc[i,0], rotation=90, color=color[i], fontsize=8, verticalalignment='top', horizontalalignment='center')

    # labels
    ax.set_xlabel("RT")
    ax.set_ylabel("intensity")
    # add lines and labels
    for i, line in enumerate(lines):
        ax.axvline(line, color="red", linestyle="--", linewidth=0.5)
        ax.text(line, 0, labels[i], rotation=90, color='red', fontsize=8, verticalalignment='top', horizontalalignment='center')

    # legend
    legend_names = [[Path(file).name, Path(file).name] for file in all_files_RT_axis_vs_mass_sum_halocarbs]
    # flatten
    legend_names = [item for sublist in legend_names for item in sublist]

    ax.legend(legend_names, loc='upper right')





    sum_over_mass_halocarbs_eval = sum_over_mass_halocarbs_eval.reshape((sum_over_mass_halocarbs_eval.shape[0], 1))
    # loop over all files and sum up
    for file in all_files_RT_axis_vs_mass_sum_halocarbs[1:]:
        # load data
        RT_axis = np.loadtxt(file, delimiter="\t", usecols=0) #x
        sum_over_mass_halocarbs = np.loadtxt(file, delimiter="\t", usecols=2) #y
        sum_over_mass_halocarbs = sum_over_mass_halocarbs.reshape((sum_over_mass_halocarbs.shape[0], 1))
        # add to sum
        RT_axis_eval = np.hstack((RT_axis_eval, RT_axis))
        sum_over_mass_halocarbs_eval = np.hstack((sum_over_mass_halocarbs_eval, sum_over_mass_halocarbs))
    # calculate mean and std dev
    sum_over_mass_halocarbs_mean = np.mean(sum_over_mass_halocarbs_eval, axis=1)
    sum_over_mass_halocarbs_std = np.std(sum_over_mass_halocarbs_eval, axis=1)
    # plot mean and std dev
    fig, ax = plt.subplots(1,1)
    ax.plot(RT_axis_eval, sum_over_mass_halocarbs_mean)
    ax.fill_between(RT_axis_eval, sum_over_mass_halocarbs_mean - sum_over_mass_halocarbs_std, sum_over_mass_halocarbs_mean + sum_over_mass_halocarbs_std, alpha=0.5)
    # add calibrant lines and labels




    # plot mass axis vs sum over RT
    fig, ax = plt.subplots(1,1)
    ax.plot(RT_axis, sum_over_mass)
    ax.plot(RT_axis, sum_over_mass_halocarbs)
    ax.plot(RT_axis, sum_over_mass_halocarbs.sum(axis=1))
    # labels
    ax.set_xlabel("RT")
    ax.set_ylabel("intensity")
    # add calibrant lines and labels
    #lines until mz = 150
    # legend
    ax.legend(["all masses", '20','34','45','46','47','59','60', "sum over halocarbons"])
        # set limit from 20 to 350
                
    calculate_total_mass_spectra = False
    RT_center_where_substance_was_found = [1898, 1839, 1852, 1995, 2005, 1978, 1881, 1898, 1903, 1538]
    all_files_where_substance_was_found = ['C:\\Users\\kaho\\Desktop\\data\\data_Empa\\results_campaign2023\\BASF_real_duplicate_2\\230311.0125.tank.1.h5',
                                            'C:\\Users\\kaho\\Desktop\\data\\data_Empa\\results_campaign2023\\BASF_single_1\\230313.1808.tank.3.h5',
                                            'C:\\Users\\kaho\\Desktop\\data\\data_Empa\\results_campaign2023\\Chemours_real_duplicate_2\\230309.2224.tank.1.h5',
                                            'C:\\Users\\kaho\\Desktop\\data\\data_Empa\\results_campaign2023\\Givaudon_real_duplicate_2_first_run_discrete_fitting\\230310.0724.tank.5.h5',
                                            'C:\\Users\\kaho\\Desktop\\data\\data_Empa\\results_campaign2023\\Givaudon_single_1\\230310.1925.tank.11.h5',
                                            'C:\\Users\\kaho\\Desktop\\data\\data_Empa\\results_campaign2023\\Schweizerhalle_real_duplicate_2\\230310.1325.tank.9.h5',
                                            'C:\\Users\\kaho\\Desktop\\data\\data_Empa\\results_campaign2023\\Schweizerhalle_single_1\\230311.1625.tank.9.h5',
                                            'C:\\Users\\kaho\\Desktop\\data\\data_Empa\\results_campaign2023\\Solvay_real_duplicate_2\\230311.0555.tank.3.h5',
                                            'C:\\Users\\kaho\\Desktop\\data\\data_Empa\\results_campaign2023\\Solvay_single_2\\230311.1025.tank.5.h5']
    

    all_files_where_substance_was_found = ['C:\\Users\\kaho\\Desktop\\data\\data_Empa\\results_campaign2023\\Chemours_real_duplicate_2\\230309.2224.tank.1.h5', 'C:\\Users\\kaho\Desktop\\data\\data_Empa\\Campaign202303\\230308.1924.tank.1.h5']
    RT_center_where_substance_was_found = [1394,1475,1506, 1565, 1897,1943,2032,2149,2861]


    calculate_total_mass_spectra = False
    if calculate_total_mass_spectra:

        # calcualate parameters 
        for i, file in enumerate(all_files_where_substance_was_found):
            # interesting RT ranges
            RT_ranges_eval = []
            for ii, RT in enumerate(RT_center_where_substance_was_found):
                RT_range_where_substance_was_found_i = [RT + RT_rel_range_where_substance_was_found[0], RT + RT_rel_range_where_substance_was_found[1]]
                RT_ranges_eval.append(RT_range_where_substance_was_found_i)

                for spectrum_type in ['EI', 'CI']:
                    data_ionimode, reader = load_data(file, ionisation_dict, spectrum_type)
                    # get mass axis
                    mass_axis = f_tofidx_to_mass(np.arange(data_ionimode.shape[1]), reader.mass_calibration_data[spectrum_type].get_parameters_for_rt(RT))
                    # get sum over RT +- RT_rel_range_where_substance_was_found
                    sum_over_RT = get_sum(data_ionimode, sum_type = 'mass', tol=0.1, range_starts=[(RT + RT_rel_range_where_substance_was_found[0])/s_per_bin], range_ends=[(RT + RT_rel_range_where_substance_was_found[1])/s_per_bin])
                    # plot mass axis vs sum over RT
                    fig, ax = plt.subplots(1,1)
                    ax.plot(sum_over_RT)
                    ax.plot(mass_axis, sum_over_RT)
                    # set limit from 20 to 350
                    ax.set_xlim([20, 350])
                    # add calibrant lines and labels
                    #for iii, line in enumerate(lines):
                    #    ax.axvline(line, color="red", linestyle="--", linewidth=0.5)
                    #    ax.text(line, 0, label[iii], rotation=90, color='red', fontsize=8, verticalalignment='top', horizontalalignment='center')
                    # save as python plot
                    filename_str = folder_to_save_intermediate_results + r"\mass_axis_vs_RT_sum_{}_{}_{}.png".format(Path(file).stem, spectrum_type, RT)
                    plt.savefig(filename_str, dpi=300, bbox_inches='tight')
                    plt.close()

                    # save data: first column mass axis, second column and following: RT
                    data_to_save_RT = data_ionimode[RT_range_where_substance_was_found_i[0]:RT_range_where_substance_was_found_i[1], : ]
                    np.savetxt(folder_to_save_intermediate_results + r"\mass_axis_vs_RT_sum_{}_{}_{}.txt".format(Path(file).stem, spectrum_type, RT), data_to_save_RT, delimiter="\t")
                    # print dimensions of mass axis vs RT sum
                    # bring mass axis and sum over RT in the same shape
                    mass_axis = mass_axis.reshape((mass_axis.shape[0], 1))
                    data_to_save_mass = np.hstack((mass_axis, sum_over_RT))
                    np.savetxt(folder_to_save_intermediate_results + r"\mass_axis_vs_RT_sum_{}_{}_{}_mass.txt".format(Path(file).stem, spectrum_type, RT), data_to_save_mass, delimiter="\t")

    stop()
    # reload x values, the mass axis and y array values, the sum over interpolated mass axis

    #sum_over_RT = get_sum(data_ionimode, sum_type = 'RT', tol=0.1, range_starts=[(RT + RT_rel_range_where_substance_was_found[0])/s_per_bin], range_ends=[(RT + RT_rel_range_where_substance_was_found[1])/s_per_bin])

    # get all files in folder
    all_files = [str(file) for file in Path(folder_to_save_intermediate_results).glob('**/*.txt')]
    # get all files that contain the mass axis and the sum over RT
    all_files_mass_axis_vs_RT_sum = [file for file in all_files if "mass_axis_vs_RT_sum" in file and "_mass.txt" not in file]
    # get all files that contain the mass axis and the sum over RT
    all_files_mass_axis_vs_RT_sum_mass = [file for file in all_files if "_mass.txt" in file]

    RT = 1839
    RT =  1565
    RT_rel_range_where_substance_was_found
    range_starts=[(RT + RT_rel_range_where_substance_was_found[0])/s_per_bin]
    range_ends=[(RT + RT_rel_range_where_substance_was_found[1])/s_per_bin]
    # get all files that are specific for one RT
    all_files_specific_RT = [file for file in all_files_mass_axis_vs_RT_sum if str(RT) in file]
    all_files_specific_RT_mass = [file for file in all_files_mass_axis_vs_RT_sum_mass if str(RT) in file]

    # prepare empty arrays for summing
    # x is mass axis of first file
    mz_gridd_eval = np.loadtxt(all_files_specific_RT_mass[0], delimiter="\t", usecols=0) #x
    RT_gridd_eval = np.arange(len(np.loadtxt(all_files_mass_axis_vs_RT_sum[0], delimiter="\t", usecols=0))) #y

    interpolated_values_list = []

    spectra_mode= "EI"
    for file in all_files_specific_RT :
        if spectra_mode in file:
            intensity = np.loadtxt(file, delimiter="\t")
            # shape of f_values is (x, y)
            RT_dim, mz_dim = np.shape(intensity) # all masses, RT indici in range
            # get mass axis
            file_mass = file.replace(".txt", "_mass.txt")
            mz_gridd = np.loadtxt(file_mass, delimiter="\t", usecols=0) #x
            # print shape
            np.shape(mz_gridd) # all masses
            RT_gridd = np.arange(int(range_starts[0]*s_per_bin),int(range_ends[0]*s_per_bin)) #y
            # Initialize an empty array to store the interpolated values
            interpolated_values = np.empty((len(RT_gridd), len(mz_gridd_eval)))

            # Loop over the `RT` dimension and perform interpolation for each `RT`
            for i in range(RT_dim):
                # Create an interpolation function for the `mz` dimension at the fixed `RT`
                interp_function = interp1d(mz_gridd, intensity[i, :], kind='linear', fill_value=0.0, bounds_error=False)

                # Interpolate the data for the given `mz_eval`
                interpolated_values[i, :] = interp_function(mz_gridd_eval)
            
            interpolated_values_list.append(interpolated_values)









    # calculate mean and std dev
    interpolated_values_list = np.array(interpolated_values_list)
    interpolated_values_mean = np.mean(interpolated_values_list, axis=0)
    # print shape
    np.shape(interpolated_values_mean) # all masses
    interpolated_values_std = np.std(interpolated_values_list, axis=0)

    # plot sum over all masses
    fig, ax = plt.subplots(1,1)
    RT_gridd_eval_plot = (RT_gridd_eval -20+RT/s_per_bin)*s_per_bin
    ax.plot(RT_gridd_eval_plot, interpolated_values_mean.sum(axis=1))
    ax.fill_between(RT_gridd_eval_plot, interpolated_values_mean.sum(axis=1) - interpolated_values_std.sum(axis=1), interpolated_values_mean.sum(axis=1) + interpolated_values_std.sum(axis=1), alpha=0.5)
    # add smoothing
    from scipy.signal import savgol_filter
    sum_over_all_masses_smoothed = savgol_filter(interpolated_values_mean.sum(axis=1), 11, 3)
    ax.plot(RT_gridd_eval_plot, sum_over_all_masses_smoothed)


    # # same for mass 
    # fig, ax = plt.subplots(1,1)
    # mz_gridd_eval = mz_gridd_eval
    # ax.plot(mz_gridd_eval, interpolated_values_mean.sum(axis=0))
    # ax.fill_between(mz_gridd_eval, interpolated_values_mean.sum(axis=0) - interpolated_values_std.sum(axis=0), interpolated_values_mean.sum(axis=0) + interpolated_values_std.sum(axis=0), alpha=0.5)
    # # add smoothing
    # from scipy.signal import savgol_filter
    # #sum_over_all_masses_smoothed = savgol_filter(interpolated_values_mean.sum(axis=0), 11, 3)
    # ax.plot(mz_gridd_eval, sum_over_all_masses_smoothed)

    # plot first third, center third and last third
    interpolated_values_mean_first_third = interpolated_values_mean[:int(np.shape(interpolated_values_mean)[0]/3),:].sum(axis=0)
    interpolated_values_mean_center_third = interpolated_values_mean[ int(np.shape(interpolated_values_mean)[0]/3):int(2*np.shape(interpolated_values_mean)[0]/3),:].sum(axis=0)
    interpolated_values_mean_last_third = interpolated_values_mean[int(2*np.shape(interpolated_values_mean)[0]/3):,:].sum(axis=0)
    interpolated_values_std_first_third = interpolated_values_std[:int(np.shape(interpolated_values_mean)[0]/3),:].sum(axis=0)
    interpolated_values_std_center_third = interpolated_values_std[int(np.shape(interpolated_values_mean)[0]/3):int(2*np.shape(interpolated_values_mean)[0]/3),:].sum(axis=0)
    interpolated_values_std_last_third = interpolated_values_std[int(2*np.shape(interpolated_values_mean)[0]/3):,:].sum(axis=0)
    fig, ax = plt.subplots(1,1)
    ax.plot(mz_gridd_eval, interpolated_values_mean_first_third)
    ax.fill_between(mz_gridd_eval, interpolated_values_mean_first_third - interpolated_values_std_first_third, interpolated_values_mean_first_third + interpolated_values_std_first_third, alpha=0.2)
    ax.plot(mz_gridd_eval, interpolated_values_mean_center_third)
    ax.fill_between(mz_gridd_eval, interpolated_values_mean_center_third - interpolated_values_std_center_third, interpolated_values_mean_center_third + interpolated_values_std_center_third, alpha=0.2)
    ax.plot(mz_gridd_eval, interpolated_values_mean_last_third)
    ax.fill_between(mz_gridd_eval, interpolated_values_mean_last_third - interpolated_values_std_last_third, interpolated_values_mean_last_third + interpolated_values_std_last_third, alpha=0.2)
    # labels
    ax.set_xlabel("mass")
    ax.set_ylabel("intensity")
    # add legend
    ax.legend(["first third", "first third uncertainty", "center third", "center third uncertainty", "last third", "last third uncertainty"])


    # same, but plot difference center - first third and  center - second third
    interpolated_values_mean_center_minus_first_third = interpolated_values_mean_center_third - interpolated_values_mean_first_third
    interpolated_values_mean_center_minus_last_third = interpolated_values_mean_center_third - interpolated_values_mean_last_third
    interpolated_values_std_center_minus_first_third = np.sqrt(interpolated_values_std_center_third**2 + interpolated_values_std_first_third**2)
    interpolated_values_std_center_minus_last_third = np.sqrt(interpolated_values_std_center_third**2 + interpolated_values_std_last_third**2)
    fig, ax = plt.subplots(1,1)
    ax.plot(mz_gridd_eval, interpolated_values_mean_center_minus_first_third)
    ax.fill_between(mz_gridd_eval, interpolated_values_mean_center_minus_first_third - interpolated_values_std_center_minus_first_third, interpolated_values_mean_center_minus_first_third + interpolated_values_std_center_minus_first_third, alpha=0.2)
    ax.plot(mz_gridd_eval, interpolated_values_mean_center_minus_last_third)
    ax.fill_between(mz_gridd_eval, interpolated_values_mean_center_minus_last_third - interpolated_values_std_center_minus_last_third, interpolated_values_mean_center_minus_last_third + interpolated_values_std_center_minus_last_third, alpha=0.2)
    # labels
    ax.set_xlabel("mass")
    ax.set_ylabel("intensity")
    # add legend
    ax.legend(["center - first third", "center - first third uncertainty", "center - last third", "center - last third uncertainty"])




    # add smoothing




    # add calibrant lines and labels
    #lines until mz = 150


    # plot mean and std dev
    from mpl_toolkits import mplot3d
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # plot 3d data
    # print shapes
    np.shape(RT_gridd_eval)
    np.shape(mz_gridd_eval)
    np.shape(interpolated_values_mean)
    surface = ax.plot_surface(RT_gridd_eval[0], mz_gridd_eval[1], interpolated_values_mean, cmap='viridis', edgecolor='none')

    # Add labels and a color bar
    ax.set_xlabel('Retention Time (RT)')
    ax.set_ylabel('Mass-to-Charge Ratio (m/z)')
    ax.set_zlabel('Intensity')
    fig.colorbar(surface, label='Intensity')

    # Show the plot
    plt.show()
    ax.plot(mz_gridd_eval, interpolated_values_mean)
    ax.fill_between(mz_gridd_eval, interpolated_values_mean - interpolated_values_std, interpolated_values_mean + interpolated_values_std, alpha=0.5)
    # add calibrant lines and labels
    #lines until mz = 150
    index_lines_below_150 = np.where(lines < 150)[0]
    for i, line in enumerate(lines[index_lines_below_150]):
        ax.axvline(line, color="red", linestyle="--", linewidth=0.5)
        ax.text(line, 0, label[i], rotation=90, color='red', fontsize=8, verticalalignment='top', horizontalalignment='center')



    interpolated_values_list_first_half = []
    interpolated_values_list_second_half = []
    interpolated_values_list_center = []

    # sum over RT
    masses_to_ignore = [185, 187]
    # get index where mz_gridd is > 23 and < 120
    mz_gridd_eval_ind = np.where((mz_gridd > 23) & (mz_gridd < 120))[0]
    # exclude also values between 56.5 and 57.5
    mz_gridd_eval_ind = np.where((mz_gridd > 23) & (mz_gridd < 120) & (mz_gridd < 56.5) | (mz_gridd > 57.5))[0]
    sum_over_all_masses = f_values[:,mz_gridd_eval_ind].sum(axis=1)
    # plot against RT_gridd
    fig, ax = plt.subplots(1,1)
    ax.plot(RT_gridd, sum_over_all_masses)
    # plot smoothed
    from scipy.signal import savgol_filter
    sum_over_all_masses_smoothed = savgol_filter(sum_over_all_masses, 11, 3)
    ax.plot(RT_gridd, sum_over_all_masses_smoothed)

    # sum over axis
    center_five_ind = np.where((RT_gridd > 1898 -2) & (RT_gridd < 1898 + 3))[0]
    left_five_ind = np.where((RT_gridd > 1898 - 7) & (RT_gridd < 1898 - 2))[0]
    right_five_ind = np.where((RT_gridd > 1898 + 3) & (RT_gridd < 1898 + 8))[0]
    f_values_sum_center = f_values[center_five_ind, :].sum(axis=0)
    f_values_sum_first_half = f_values[left_five_ind,:].sum(axis=0)
    f_values_sum_second_half = f_values[right_five_ind,:].sum(axis=0)
    # print mass axis and summed up values
    fig, ax = plt.subplots(1,1)
    ax.plot(mz_gridd, f_values_sum_first_half)
    ax.plot(mz_gridd, f_values_sum_second_half)
    ax.plot(mz_gridd, f_values_sum_center)
    # add labels
    ax.set_xlabel("mass")
    ax.set_ylabel("intensity")
    # add legend
    ax.legend(["first half", "second half", "center"])




    # calculate mean and std dev
    interpolated_values_list = np.array(interpolated_values_list)
    interpolated_values_mean = np.mean(interpolated_values_list, axis=0)
    interpolated_values_std = np.std(interpolated_values_list, axis=0)

    # plot mean and std dev
    fig, ax = plt.subplots(1,1)
    ax.plot(mz_gridd_eval, interpolated_values_mean)
    ax.fill_between(mz_gridd_eval, interpolated_values_mean - interpolated_values_std, interpolated_values_mean + interpolated_values_std, alpha=0.5)
    # add calibrant lines and labels
        

    # interpolated_values now contains the interpolated values at the arbitrary coordinates
    print(interpolated_values)







    # load data
    data_EI, reader_EI = load_data(file, ionisation_dict, 'EI')
    reader_EI.mass_calibration_data['EI'].get_parameters_for_rt_nearest_neighbours() #updated myguversion was necessary to get the correct parameters
    reader_EI.mass_calibration_data['EI'].get_parameters_for_rt(1455) #updated myguversion was necessary to get the correct parameters
    mass_axis_EI = reader_EI.get_mass_axis_of(1898/0.2)
    # get parameters for mass axis
    calib_params_EI = reader_EI.mass_calibration_data['EI'].get_parameters_for_rt(1898)
    # same for CI
    data_CI, reader_CI = load_data(file, ionisation_dict, 'CI')
    #get calibration parameters from mc file
    #mass_cal_param_list_for_RT = reader_EI.mass_calibration_data['EI'].mass_cal_parameters
    # get calibration parameters from hdf5 file


    # def calibrate_via_surface_fit(mass_cal_param_list_for_RT):
    #     data = np.array(mass_cal_param_list_for_RT)
    #     # convert to numpy array
    #     RT = np.array([data[i][0] for i in range(len(data))])
    #     p1 = np.array([data[i][1] for i in range(len(data))])
    #     p2 = np.array([data[i][2] for i in range(len(data))])
    #     p3 = np.array([data[i][3] for i in range(len(data))])

    #     from sklearn.preprocessing import PolynomialFeatures
    #     from sklearn.linear_model import LinearRegression

    #     # Sample calibration data (replace with your actual data)
    #     times = RT
    #     params = np.array([[a, b, c] for a, b, c in zip(p1, p2, p3)])

    #     # Fit a 3D polynomial regression model
    #     poly = PolynomialFeatures(degree=1)  # should be a simple linear regression
    #     X = poly.fit_transform(times[:, np.newaxis])
    #     model = LinearRegression()
    #     model.fit(X, params)

    #     return model, poly
    
    # def get_calibration_paramters(model, poly, RT_eval):
    #     # Evaluate the model at time 't'
    #     X_new = poly.transform(np.array([[RT_eval]]))
    #     estimated_params = model.predict(X_new)

    #     print("Estimated Parameters at t =", RT_eval, ":", estimated_params)

    #     return estimated_params
    
    #model, poly = calibrate_via_surface_fit(mass_cal_param_list_for_RT)
    #get_calibration_paramters(model, poly, 1455)


    #save mass_axis_EI to file and intensities for 20 datapoints around 1898 to txtfile
  


    # get sum over RT +- 2s
    sum_over_RT_EI = get_sum(data_EI, sum_type = 'mass', tol=0.1, range_starts=[(1898-0.5)/0.2], range_ends=[(1898+0.5)/0.2])

    # same for CI
    mass_axis_CI = reader_CI.get_mass_axis_of(1898/0.2)
    # get parameters for mass axis
    calib_params_CI = reader_CI.mass_calibration_data['CI'].get_parameters_for_rt(1898)

    # get sum over RT +- 2s
    sum_over_RT_CI = get_sum(data_CI, sum_type = 'mass', tol=0.1, range_starts=[(1898-0.5)/0.2], range_ends=[(1898+0.5)/0.2])


    # plot mass axis vs sum over RT
    fig, ax = plt.subplots(1,1)
    ax.plot(mass_axis_EI, sum_over_RT_EI)
    # add CI
    ax.plot(mass_axis_CI, sum_over_RT_CI)
    # add calibrant lines
    calibrant_data_path = r"c:\Users\kaho\Desktop\data\data_Empa\combinatorial_data\formulas_CF_radicals_True.txt"
    calibrant_data = pd.read_csv(calibrant_data_path, sep="\t", header=None, names=["formula", "mass"])
    lines = calibrant_data['mass']
    label = calibrant_data['formula']
    # add calibrant lines and labels
    for i, line in enumerate(lines):
        ax.axvline(line, color="red", linestyle="--", linewidth=0.5)
        ax.text(line, 0, label[i], rotation=90, color='red', fontsize=8, verticalalignment='top', horizontalalignment='center')


    #add metadata as title
    plt.title("RT range: {} - {}".format(1897, 1899))

    # plot largest 50 masses over RT
    fig, ax = plt.subplots(1,1)
    ax.plot(data_EI[:, int(1800/0.2):int(1950/0.2)].sum(axis=1))
    #add metadata as title
    plt.title("RT range: {} - {}".format(1897, 1899))






    #mass_cal_by_RT_index = reader_EI.get_mass # CAREFUL!: mygu version only uses left and right neighbour
    dir(reader_EI.mass_calibration_data['EI'])


    sumd = get_sum(data_EI, sum_type = 'mass', tol=0.1)

    # PLOT RT SUM


    # sum all rt_sum data
    #fig, ax = plt.subplots(1,1)
    #rt_sum_all = rt_sum.sum(axis=1)
    # RT = index*s_per_bin + offset
    #index =  np.arange(rt_sum.shape[0])  
    #s_per_bin = float(1/5) # 5Hz
    #plt.plot(index*s_per_bin, np.array(rt_sum))


    # PLOT HIGHEST MASS PEAKS

    range_starts, range_ends = get_ranges_for_summing(data_EI[:,1890*0.2], sum_type="RT", threshold_median_multiples=2, thresh_min_range_len= 3)
    print(np.array(range_ends)-np.array(range_starts))
    # eliminiate ranges that are too small
    # print mass spectra for each range
    fig, ax = plt.subplots(1,1)
    range_center = 0.5*(range_starts + range_ends)
    line_type = ["-", "--", "-.", ":"]*100
    s_per_bin = float(1/5) # 5Hz


    for i, range_start in enumerate(range_starts):
        range_end = range_ends[i]
        rt_sum = get_sum(data_EI, sum_type="RT", range_starts=[range_start], range_ends=[range_end], tol=0.1)
        #plt.plot(rt_sum, label="RT range: "+str(range_center[i]))
        rt_sum_x = np.arange(rt_sum.shape[0])*s_per_bin
        plt.plot(rt_sum_x , rt_sum, label=round(f_tofidx_to_mass(range_center[i], reader_EI.mass_calibration_data['EI'].get_parameters_for_rt(1898))), linestyle = line_type[i])
    # add total sum
    rt_sum_all = get_sum(data_EI, sum_type="RT")
    rt_sum_all_x = np.arange(rt_sum_all.shape[0])*s_per_bin
    plt.plot(rt_sum_all_x, rt_sum_all, label="total sum")
    # plot horizontal labels
    # if color palette is repeating change the line type
    plt.legend(bbox_to_anchor=(0.05, 1), loc='upper left', borderaxespad=0.0, fontsize=6, ncol=3)
    plt.show()

    # ONE COULD IMPLEMENT A RT TIME SEPARATION HERE


    # draw mass plot with ranges
    fig, ax = plt.subplots(1,1)
    plt.plot(data.sum(axis=0))
    for i, range_start in enumerate(range_starts):
        range_end = range_ends[i]
        plt.axvline(range_start, color="red", linestyle="--", linewidth=0.5)
        plt.axvline(range_end, color="lightgreen", linestyle="--", linewidth=0.5)
        # plot horizontal line at zero
        plt.axhline(0, color="black", linewidth=0.5)

  

    print(np.array(range_ends)-np.array(range_starts))

    
    range_start, range_end = get_ranges_for_summing(data, sum_type="RT", threshold_median_multiples=5000, thresh_min_range_len= 3)
    rt_sum = get_sum(data, sum_type="RT", range_starts=[range_start], range_ends=[range_end], tol=0.1)
    get_calibrant_peaks(data, start_or_end="end", factor_median_for_threshold=10, plot=True)
    get_calibrant_peaks(data, start_or_end="start", factor_median_for_threshold=10, plot=True)


    fig, ax = plt.subplots(1,1)
    plt.plot(data.sum(axis=0)) # axis 0: mass, 1: RT
    fig, ax = plt.subplots(1,1)
    plt.plot(rt_sum)




    stop()

    data_CI, reader_CI = load_data(file, ionisation_dict, 'CI')

    
    # get mass cailbration from hdf5 file
    # manual read out of params
    # Myriams calibration
    #reader.mass_calibration_data['EI'].mass_cal_parameters
    # TOFWerk calibration
    #reader.mass_calibration_data['EI'].get_parameters_for_rt(125)
    #reader.mass_cal_tofwerk
    #reader.rt_calibration_data
    


    # get mass calibration from hdf5 file
    #read_mass_cal_par(reader)
    # i = p1*np.sqrt(mz) + p2
    p1 = 2835.23579263574
    p2 = -4049.42781522328

    #RT calibration
    s_per_bin = float(1/5) # 5Hz
    offset = 0 

    

    data = data_EI
    ranges_start_example = [600]
    ranges_end_example = [700]
    mass_spec_df= get_mass_spec_for_RT_range(data, ranges_start = ranges_start_example, ranges_end = ranges_end_example, rt_index_to_rt_mode = 0, rt_index_to_rt_params=[float(1/5),0], mass_index_to_mass_mode = 0, mass_index_to_mass_params=[2835.23579263574, -4049.42781522328], range_tol=0.1, range_mode='abs', sum_x='mass', normalize=False) 

    fig, ax = plt.subplots(1,1)
    ax.plot(mass_spec_df['mass'], mass_spec_df['intensity'])
    #add metadata as title
    plt.title("RT range: {} - {}".format(ranges_start_example, ranges_end_example   ))


    ranges_start_example = [20]
    ranges_end_example = [150]
    RT_spec_df = get_RT_spec_for_mass_range(data, ranges_start = None, ranges_end = None, rt_index_to_rt_mode = 0, rt_index_to_rt_params=[float(1/5),0], mass_index_to_mass_mode = 0, mass_index_to_mass_params=[2835.23579263574, -4049.42781522328], range_tol=0.1, range_mode='abs', sum_x='RT', normalize=False)
    fig, ax = plt.subplots(1,1)
    ax.plot(RT_spec_df['RT'], RT_spec_df['intensity'])
    #add metadata as title
    plt.title("mass range: {} - {}".format(ranges_start_example, ranges_end_example))


    # add compounds vertically written to the plot
    #read_determined_peaks_data
    path_det_peak_data = Path(r"c:\Users\kaho\polybox\Presentations\PythonPlots\peaks_det_230211.0125.tank.1.xlsx")
    det_peak_data = pd.read_excel(path_det_peak_data)
    # delete all rows which are Nan in first column
    det_peak_data = det_peak_data[det_peak_data.iloc[:,0].notna()]
    # delete rows where first column is "no usable data"
    det_peak_data = det_peak_data[det_peak_data.iloc[:,0] != "no usable data"]
    # split second column by ;\n in two columns
    det_peak_data[['RT_s', 'RT_min']] = det_peak_data.iloc[:,1].str.split(';\n',  expand=True)
    pos_label = det_peak_data['RT_s']
    # convert to float
    pos_label = pos_label.astype(float)
    label_name = det_peak_data['eval_status'].astype(str)
    # convert to np.array
    pos_label = np.array(pos_label)
    label_name = np.array(label_name)
    # add to plot as vertical labels at 0
    for i, label in enumerate(label_name):
        ax.text(pos_label[i], 0, label, rotation=90, color='red', fontsize=8, verticalalignment='top', horizontalalignment='center')

    # get maximum intensity




    # get first 4 and last 4 maxima (spikes)
    # this can be extracted from TPSscripting data and comparison to the peaks found in peakfinder
    # the calibrant trigger is at:
    calibrant_theo = [120,125,130,135,3065,3070,3075,3080]
    # number of calibrant peaks can be extracted as well

    p1 = 2835.23579263574
    p2 = -4049.42781522328


    # for fit:
    # Initial guess for p1 and p2
    p1_init = 2835.23579263574
    p2_init = -4049.42781522328
    x0 = np.array([p1_init, p2_init])
    # Define parameter bounds
    p1_bounds = (p1_init - 0.002*p1_init, p1_init + 0.002*p1_init)
    p2_bounds = (p2_init - 0.002*p2_init, p2_init + 0.002*p2_init)
    # sort bounds
    if p1_bounds[0] > p1_bounds[1]:
        p1_bounds = (p1_bounds[1], p1_bounds[0])
    if p2_bounds[0] > p2_bounds[1]:
        p2_bounds = (p2_bounds[1], p2_bounds[0])
    bounds = [p1_bounds, p2_bounds]

    # Define initial guess
    params_initial = np.array([(p1_bounds[0] + p1_bounds[1])/2, (p2_bounds[0] + p2_bounds[1])/2])



    data_CI_sum = get_sum(data_CI, "RT")
    # plot vs mass calibrated data

    
    # get sum:
    data_sum = get_sum(data, "mass")
    # plot vs mass calibrated data:
    fig, ax = plt.subplots(1,1)
    plt.plot(index_to_mass(range(len(data_sum)), p1=p1, p2 = p2), data_sum)
    rt_sum = get_sum(data, sum_type="RT", range_starts=[0, 3000], range_ends=[200, 400], tol=0.1)
    mass_sum = get_sum(data, sum_type="mass", range_starts=[10, 30], range_ends=[20, 40], tol=0.1)

    # sum all rt_sum data
    rt_sum_all = rt_sum.sum(axis=1)
    plt.plot(rt_sum_all)

    #fig, ax = plt.subplots(1,1)
    #plt.plot(rt_sum.sum)

    # get ranges for RT by choosing only the highest mass peaks


    rt_sum = get_sum(data, sum_type="RT", range_starts=[0], range_ends=[22248], tol=5)
    mass_sum = get_sum(data, sum_type="mass", tol=0.1)

    fig, ax = plt.subplots(1,1)
    plt.plot(mass_sum)

    # PLOT RT SUM


    # sum all rt_sum data
    fig, ax = plt.subplots(1,1)
    #rt_sum_all = rt_sum.sum(axis=1)
    # RT = index*s_per_bin + offset
    index =  np.arange(rt_sum.shape[0])  
    plt.plot(index*s_per_bin, np.array(rt_sum))


    # PLOT HIGHEST MASS PEAKS

    range_starts, range_ends = get_ranges_for_summing(data, sum_type="RT", threshold_median_multiples=2, thresh_min_range_len= 3)
    print(np.array(range_ends)-np.array(range_starts))
    # eliminiate ranges that are too small
    # print mass spectra for each range
    fig, ax = plt.subplots(1,1)
    range_center = 0.5*(range_starts + range_ends)
    line_type = ["-", "--", "-.", ":"]*100


    for i, range_start in enumerate(range_starts):
        range_end = range_ends[i]
        rt_sum = get_sum(data, sum_type="RT", range_starts=[range_start], range_ends=[range_end], tol=0.1)
        #plt.plot(rt_sum, label="RT range: "+str(range_center[i]))
        plt.plot(rt_sum, label=round(index_to_mass(range_center[i], p1=p1, p2=p2)), linestyle = line_type[i])
    # plot horizontal labels
    # if color palette is repeating change the line type
    plt.legend(bbox_to_anchor=(0.05, 1), loc='upper left', borderaxespad=0.0, fontsize=6, ncol=3)
    plt.show()

    # ONE COULD IMPLEMENT A RT TIME SEPARATION HERE


    # draw mass plot with ranges
    fig, ax = plt.subplots(1,1)
    plt.plot(data.sum(axis=0))
    for i, range_start in enumerate(range_starts):
        range_end = range_ends[i]
        plt.axvline(range_start, color="red", linestyle="--", linewidth=0.5)
        plt.axvline(range_end, color="lightgreen", linestyle="--", linewidth=0.5)
        # plot horizontal line at zero
        plt.axhline(0, color="black", linewidth=0.5)

  

    print(np.array(range_ends)-np.array(range_starts))

    
    range_start, range_end = get_ranges_for_summing(data, sum_type="RT", threshold_median_multiples=5000, thresh_min_range_len= 3)
    rt_sum = get_sum(data, sum_type="RT", range_starts=[range_start], range_ends=[range_end], tol=0.1)
    get_calibrant_peaks(data, start_or_end="end", factor_median_for_threshold=10, plot=True)
    get_calibrant_peaks(data, start_or_end="start", factor_median_for_threshold=10, plot=True)


    fig, ax = plt.subplots(1,1)
    plt.plot(data.sum(axis=0)) # axis 0: mass, 1: RT
    fig, ax = plt.subplots(1,1)
    plt.plot(rt_sum)