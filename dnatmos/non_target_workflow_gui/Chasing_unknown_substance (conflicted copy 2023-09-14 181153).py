




# excecute if called as main
from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from data_analysis_utils import formula_by_name, get_calibrant_peaks, get_mass_spec_for_RT_range, get_ranges_for_summing, get_sum, index_to_mass, load_data, get_RT_spec_for_mass_range
from alpinac.utils_data_extraction import f_tofidx_to_mass
from scipy.interpolate import griddata
from scipy.interpolate import interp1d


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
    all_files_where_substance_was_found = all_files[:-3]
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

    recalculate = False
    if recalculate:
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

    # get all files in folder
    all_files = [str(file) for file in Path(folder_to_save_intermediate_results).glob('**/*.txt')]
    # get all files that contain the mass axis and the sum over RT
    all_files_mass_axis_vs_RT_sum = [file for file in all_files if "mass_axis_vs_RT_sum" in file and "_mass.txt" not in file]
    # get all files that contain the mass axis and the sum over RT
    all_files_mass_axis_vs_RT_sum_mass = [file for file in all_files if "_mass.txt" in file]

    RT = 1898
    RT_rel_range_where_substance_was_found
    range_starts=[(RT + RT_rel_range_where_substance_was_found[0])/s_per_bin]
    range_ends=[(RT + RT_rel_range_where_substance_was_found[1])/s_per_bin]
    # get all files that are specific for one RT
    all_files_specific_RT = [file for file in all_files_mass_axis_vs_RT_sum if str(RT) in file]
    all_files_specific_RT_mass = [file for file in all_files_mass_axis_vs_RT_sum_mass if str(RT) in file]

    # prepare empty arrays for summing
    # x is mass axis of first file
    mz_gridd_eval = np.loadtxt(all_files_specific_RT[0], delimiter="\t", usecols=0) #x
    RT_gridd_eval = np.arange(len(np.loadtxt(all_files_mass_axis_vs_RT_sum[0], delimiter="\t", usecols=0))) #y
    interpolated_values_list_first_half = []
    interpolated_values_list_second_half = []
    interpolated_values_list_center = []
    for file in all_files_specific_RT :
        f_values = np.loadtxt(file, delimiter="\t")
        # shape of f_values is (x, y)
        np.shape(f_values) # all masses, RT indici in range
        # get mass axis
        file_mass = file.replace(".txt", "_mass.txt")
        mz_gridd = np.loadtxt(file_mass, delimiter="\t", usecols=0) #x
        # print shape
        np.shape(mz_gridd) # all masses
        RT_gridd = np.arange(int(range_starts[0]*s_per_bin),int(range_ends[0]*s_per_bin)) #y
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


        # Initialize an empty array to store the interpolated values
        interpolated_values = np.empty((len(RT_gridd), len(mz_gridd_eval)))

        # Loop over the `RT` dimension and perform interpolation for each `RT`
        for i in range(len(RT)):
            # Create an interpolation function for the `mz` dimension at the fixed `RT`
            interp_function = interp1d(mz, intensity[i, :], kind='linear', fill_value=0.0, bounds_error=False)

            # Interpolate the data for the given `mz_eval`
            interpolated_values[i, :] = interp_function(mz_gridd_eval)


        # interpolate
        
        # Create a meshgrid of the original grid coordinates
        RT_mesh, mz_mesh = np.meshgrid(RT_gridd, mz_gridd, indexing='ij')
        intensity = f_values

        # Flatten the original grid coordinates and intensities
        RT_flat = RT_mesh.flatten()
        mz_flat = mz_mesh.flatten()
        intensity_flat = intensity.flatten()

        # Create a meshgrid of the new grid coordinates
        RT_eval_mesh, mz_eval_mesh = np.meshgrid(RT_gridd_eval, mz_gridd_eval, indexing='ij')

        # Flatten the new grid coordinates
        RT_eval_flat = RT_eval_mesh.flatten()
        mz_eval_flat = mz_eval_mesh.flatten()

        # Perform linear interpolation using scipy.interpolate.griddata
        interpolated_values = griddata((RT_flat, mz_flat), intensity_flat,
                                        (RT_eval_flat, mz_eval_flat), method='linear', fill_value=0.0)

        # Reshape the interpolated values to match the shape of RT_eval and mz_eval
        interpolated_values = interpolated_values.reshape(RT_gridd_eval.shape[0], mz_gridd_eval.shape[0])

        # Now, interpolated_values contains the interpolated data on the new grid defined by RT_eval and mz_eval.




        #int_2D = np.loadtxt(file, delimiter="\t") #y
        # get dimensions of 2D array
        #dim_2D = np.shape(int_2D)
        #RT_gridd = np.arange(dim(np.loadtxt(file, delimiter="\t")).shape()) #y
        # generate x gridd
        x_mesh, y_mesh = np.meshgrid(mz_gridd, RT_gridd)
        # evaluate



        # Perform linear interpolation using scipy.interpolate.griddata
        interpolated_values = griddata((x_mesh.flatten(), y_mesh.flatten()), f_values_sum_center.flatten(),
                                        (mz_gridd_eval, RT_gridd_eval), method='linear', fill_value=0.0)
        interpolated_values_list.append(interpolated_values)

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