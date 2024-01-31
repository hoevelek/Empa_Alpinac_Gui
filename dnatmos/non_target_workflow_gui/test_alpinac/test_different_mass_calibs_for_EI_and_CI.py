from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from data_analysis_utils import load_data, get_sum, get_ranges_for_summing
file = Path(r"C:\Users\kaho\Desktop\data\data_Empa\Campaign202303\230401.1747.air.3.h5")
ei_seg = (3,4)
ci_seg = (1,2)
ionisation_dict = {"EI": ei_seg,"CI": ci_seg}
# These are some masses that we don't want as they have maybe a column bleed
masses_to_ignore = [185, 187]
data_EI = load_data(file, ionisation_dict, 'EI')
data_CI = load_data(file, ionisation_dict, 'CI')


# get sum:
data_sum_EI = get_sum(data_EI, "mass",range_starts=[1467], range_ends=[1467], tol=1)

# get sum:
data_sum_CI = get_sum(data_CI, "mass",range_starts=[1467], range_ends=[1467], tol=1)
# plot:
fig, ax = plt.subplots(1,1)
plt.plot(data_sum_EI[0])
plt.plot(-data_sum_CI[0])

data = data_sum_EI


rt_sum = get_sum(data, sum_type="RT", range_starts=[100, 300], range_ends=[200, 400], tol=0.1)
mass_sum = get_sum(data, sum_type="mass", range_starts=[10, 30], range_ends=[20, 40], tol=0.1)

# sum all rt_sum data
rt_sum_all = rt_sum.sum(axis=1)
plt.plot(rt_sum_all)

fig, ax = plt.subplots(1,1)
plt.plot(rt_sum.sum)

# get ranges for RT by choosing only the highest mass peaks


rt_sum = get_sum(data, sum_type="RT", range_starts=[22248], range_ends=[22248], tol=5)
mass_sum = get_sum(data, sum_type="mass", tol=0.1)

fig, ax = plt.subplots(1,1)
plt.plot(mass_sum)

# sum all rt_sum data
fig, ax = plt.subplots(1,1)
#rt_sum_all = rt_sum.sum(axis=1)
plt.plot(rt_sum)


    
range_starts, range_ends = get_ranges_for_summing(data, sum_type="RT", threshold_median_multiples=2, thresh_min_range_len= 3)
print(np.array(range_ends)-np.array(range_starts))
# eliminiate ranges that are too small
# print mass spectra for each range
fig, ax = plt.subplots(1,1)
range_center = 0.5*(range_starts + range_ends)

for i, range_start in enumerate(range_starts):
    range_end = range_ends[i]
    rt_sum = get_sum(data, sum_type="RT", range_starts=[range_start], range_ends=[range_end], tol=0.1)
    plt.plot(rt_sum, label="RT range: "+str(range_center[i]))
plt.legend()
plt.show()

# draw mass plot with ranges
fig, ax = plt.subplots(1,1)
plt.plot(data.sum(axis=0))
for i, range_start in enumerate(range_starts):
    range_end = range_ends[i]
    plt.axvline(range_start, color="red", linestyle="--", linewidth=0.5)
    plt.axvline(range_end, color="lightgreen", linestyle="--", linewidth=0.5)
    # plot horizontal line at zero
    plt.axhline(0, color="black", linewidth=0.5)