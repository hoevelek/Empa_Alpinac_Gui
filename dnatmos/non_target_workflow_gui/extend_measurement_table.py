import pandas as pd
from datetime import datetime

path_orig = r"g:\503_Themen\Klima\ecToFKampagne\Gemessene_Proben\samplelist.xlsx"
path_meas_campaign = r"g:\503_Themen\Klima\ecToFKampagne\Gemessene_Proben\Measurementplan_ecTOF_230405.xlsx"
path_extended = r"g:\503_Themen\Klima\ecToFKampagne\Gemessene_Proben\samplelist_extended.xlsx"
path_data = r"c:\Users\kaho\Desktop\data\data_Empa\Campaign202303"
# if path_tmp_processing_data does not exist create it
if not os.path.exists(path_tmp_processing_data):
    os.makedirs(path_tmp_processing_data)

# Read original file
df_orig = pd.read_excel(path_orig)
# print colnames
print(df_orig.columns)

# Read measurement campaign file
df_meas_campaign = pd.read_excel(path_meas_campaign)
# print colnames
print(df_meas_campaign.columns)


# create dictionary df_meas_campaign['Tank'] by df_meas_campaign['Filling']
dict_tank = dict(zip(df_meas_campaign['Filling'], df_meas_campaign['Tank']))

# modify df_meas_campaign['Company nearby'] such that first occurence is with ' 1', second with ' 2', ...
df_meas_campaign['Company nearby'] = df_meas_campaign['Company nearby'].astype(str) + "_" + df_meas_campaign.groupby('Company nearby').cumcount().add(1).astype(str)
# remove nan strings
df_meas_campaign['Company nearby'] = df_meas_campaign['Company nearby'].str.replace('nan', '')
# remove .0 strings
df_meas_campaign['Company nearby'] = df_meas_campaign['Company nearby'].str.replace('.0', '')


# create dictionary df_meas_campaign['Company nearby'] by df_meas_campaign['Filling']. 
dict_company_nearby = dict(zip(df_meas_campaign['Filling'], df_meas_campaign['Company nearby']))

# same for sample type
dict_sample_type = dict(zip(df_meas_campaign['Filling'], df_meas_campaign['Sample type']))

# same for site/filling
dict_site_filling = dict(zip(df_meas_campaign['Filling'], df_meas_campaign['Site / Filling']))

# extend df_orig with new columns
df_orig['Tank'] = df_orig['sample'].map(dict_tank)
df_orig['Company nearby'] = df_orig['sample'].map(dict_company_nearby)
df_orig['Sample type'] = df_orig['sample'].map(dict_sample_type)
df_orig['Site / Filling'] = df_orig['sample'].map(dict_site_filling)

# write extended df_orig to excel
df_orig.to_excel(path_extended, index=False)

# go to 'Company nearby' and get the string df_orig[date] + df_orgig[time] + df_orig['tank'] joined  by "."
# if 'Company nearby' endswith '2' add df_orig[date] + df_orgig[time] + df_orig['tank'] to dict_for_processing_2nd
# if 'Company nearby' endswith '1' add df_orig[date] + df_orgig[time] + df_orig['tank'] to dict_for_processing_1st
# if 'Company nearby' endswith '3' add df_orig[date] + df_orgig[time] + df_orig['tank'] to dict_for_processing_3rd
list_for_processing_2nd = []
list_for_processing_1st = []
list_for_processing_3rd = []
list_for_processing_4th = []
company_nearby_2nd = []
company_nearby_1st = []
company_nearby_3rd = []
company_nearby_4th = []

for index, row in df_orig.iterrows():
    time_zero_padded = str(row['time']).zfill(4)
    if str(row['Company nearby']).endswith('2'):
        #dict_for_processing_2nd[row['Company nearby']] = datetime.strptime(str(row['date']) + str(row['time']), "%y%m%d%H%M")
        print(row['Company nearby'])
        listelem = '.'.join([str(row['date']),time_zero_padded,str(row['type'])])
        list_for_processing_2nd.append(listelem)
        company_nearby_2nd.append(row['Company nearby'])
    
    elif str(row['Company nearby']).endswith('1'):
        #dict_for_processing_1st[row['Company nearby']] = datetime.strptime(str(row['date']) + str(row['time']), "%y%m%d%H%M")
        listelem = '.'.join([str(row['date']),time_zero_padded,str(row['type'])])
        list_for_processing_1st.append(listelem)
        company_nearby_1st.append(row['Company nearby'])

    elif str(row['Company nearby']).endswith('3'):
        listelem = '.'.join([str(row['date']),time_zero_padded,str(row['type'])])
        list_for_processing_3rd.append(listelem)
        company_nearby_3rd.append(row['Company nearby'])
        #dict_for_processing_3rd[row['Company nearby']] = datetime.strptime(str(row['date']) + str(row['time']), "%y%m%d%H%M")

    elif str(row['Company nearby']).endswith('4'):
        listelem = '.'.join([str(row['date']),time_zero_padded,str(row['type'])])
        list_for_processing_4th.append(listelem)
        company_nearby_4th.append(row['Company nearby'])
        #dict_for_processing_3rd[row['Company nearby']] = datetime.strptime(str(row['date']) + str(row['time']), "%y%m%d%H%M")


# join lists and companys to pd dataframe
df_2nd = pd.DataFrame({'Company nearby': company_nearby_2nd, 'date.time.type': list_for_processing_2nd})
df_1st = pd.DataFrame({'Company nearby': company_nearby_1st, 'date.time.type': list_for_processing_1st})
df_3rd = pd.DataFrame({'Company nearby': company_nearby_3rd, 'date.time.type': list_for_processing_3rd})
df_4th = pd.DataFrame({'Company nearby': company_nearby_4th, 'date.time.type': list_for_processing_4th})


# get list_for_processing_2nd value and remove the first occurence of same 'company nearby' value

# get list_for_processing_2nd value and keep only each second occurence of same 'company nearby' value
duplicates_mask = df_2nd.drop_duplicates(subset=['Company nearby'], keep='first')
second_occurence_mask = (df_2nd.groupby(df_2nd['Company nearby']).cumcount() == 1)
filtered_df = df_2nd[second_occurence_mask]
# generate a folder for processing_2nd_sample and copy list there
filtered_df.to_csv(path_tmp_processing_data + r"\list_for_processing_2nd_sample.csv", index=False)
path_tmp_processing_data = r"g:\503_Themen\Klima\ecToFKampagne\Gemessene_Proben\processing_2nd_sample"


# same for list_for_processing_1st
duplicates_mask = df_1st.drop_duplicates(subset=['Company nearby'], keep='first')
first_occurence_mask = (df_1st.groupby(df_1st['Company nearby']).cumcount() == 1)
filtered_df = df_1st[first_occurence_mask]
# generate a folder for processing_1st_sample and copy list there
filtered_df.to_csv(path_tmp_processing_data + r"\list_for_processing_1st_sample.csv", index=False)
path_tmp_processing_data = r"G:\503_Themen\Klima\kaho\preprocessed_data_campaign2023\first_sample_per_location"

# same for list_for_processing_3rd
duplicates_mask = df_3rd.drop_duplicates(subset=['Company nearby'], keep='first')
third_occurence_mask = (df_3rd.groupby(df_3rd['Company nearby']).cumcount() == 1)
filtered_df = df_3rd[third_occurence_mask]
# generate a folder for processing_3rd_sample and copy list there
filtered_df.to_csv(path_tmp_processing_data + r"\list_for_processing_3rd_sample.csv", index=False)
path_tmp_processing_data = r"G:\503_Themen\Klima\kaho\preprocessed_data_campaign2023\third_sample_per_location"

# same for list_for_processing_4th
duplicates_mask = df_4th.drop_duplicates(subset=['Company nearby'], keep='first')
fourth_occurence_mask = (df_4th.groupby(df_4th['Company nearby']).cumcount() == 1)
filtered_df = df_4th[fourth_occurence_mask]
# generate a folder for processing_4th_sample and copy list there
filtered_df.to_csv(path_tmp_processing_data + r"\list_for_processing_4th_sample.csv", index=False)
path_tmp_processing_data = r"G:\503_Themen\Klima\kaho\preprocessed_data_campaign2023\fourth_sample_per_location"






# copy h5 file with starts with filtered df[date.time.type] from folder path_data to folder sample_folder
import os
import shutil

tanknames = filtered_df['date.time.type']
#tanknames = ['230311.1625.tank']

for file_process_name in tanknames:
    print(file_process_name)    
    for file in os.listdir(path_data):
        # data
        if file.startswith(file_process_name) and file.endswith('.h5'):
            shutil.copy(path_data + "\\" + file, path_tmp_processing_data)
        # masscalibration data ending with .mc.txt
        elif file.startswith(file_process_name) and file.endswith('_mc.txt'):
            shutil.copy(path_data + "\\" + file, path_tmp_processing_data)

    # copy directory with name filtered df[date.time.type] from folder path_data to folder sample_folder 
    for file in os.listdir(path_data):
        if file.startswith(file_process_name) and os.path.isdir(path_data + "\\" + file):
            if not os.path.exists(path_tmp_processing_data + "\\" + file):
                shutil.copytree(path_data + "\\" + file, path_tmp_processing_data + "\\" + file)
        















