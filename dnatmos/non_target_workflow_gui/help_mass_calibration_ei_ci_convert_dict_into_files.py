import os

# EI part
rt_EI = {} 

#rt_EI['CH2Cl2'] = 1642.1
rt_EI['Benzene'] = 1892.4
rt_EI['Toluene'] = 2213.0
rt_EI["Trichlorofluoromethane"]= 1596.1
rt_EI["Dichlorodifluoromethane"] = 1455.7
rt_EI["Tetrachloromethane"] = 1739.1

directory = r"C:\Users\kaho\Desktop\data\data_Empa\test_mass_cal"
# different mass cal settings 
rt_EI['Benzene'] = 1908
rt_EI['Toluene'] = 2216.0
rt_EI["Trichlorofluoromethane"]= 1597
rt_EI["Dichlorodifluoromethane"] = 1468
rt_EI["Tetrachloromethane"] = 1739.1

directory = r"G:\503_Themen\Klima\ecToFKampagne\Gemessene_Proben\processing_2nd_sample"
# company nearby	date.time.type
# 2	Chemours (Dordrecht)_2	230309.2224.tank
# 4	Givaudon (Kemptal)_2	230310.0855.tank
# 7	Schweizerhalle (Basel)_2	230310.1325.tank
# 8	BASF (Ludwigshafen)_2	230311.0125.tank
# 11	Solvay (Bad Wimpfen)_2	230311.0555.tank
# 14	Eawag ARA (Dübendorf)_2	230315.0350.tank

#1 
rt_EI['Benzene'] = 1908
rt_EI['Toluene'] = 2216.0
rt_EI["Trichlorofluoromethane"]= 1597
rt_EI["Dichlorodifluoromethane"] = 1468
rt_EI["Tetrachloromethane"] = 1739.1

dict_rt_EI_for_file = {}
dict_rt_EI_for_file['230309.2224.tank'] = rt_EI

#4
rt_EI['Benzene'] = 1908.143286
rt_EI['Toluene'] = 2215.6
rt_EI["Trichlorofluoromethane"]= 1597.15
rt_EI["Dichlorodifluoromethane"] = 1467.646519
rt_EI["Tetrachloromethane"] = 1740.09606

dict_rt_EI_for_file['230310.1325.tank'] = rt_EI

#6
rt_EI['Benzene'] = 1908.479381
rt_EI['Toluene'] = 2215.94723
rt_EI["Trichlorofluoromethane"]= 1597.372393
rt_EI["Dichlorodifluoromethane"] = 1468.148516
rt_EI["Tetrachloromethane"] = 1740.769826

dict_rt_EI_for_file['230310.1325.tank'] = rt_EI

#8
rt_EI['Benzene'] = 1908.479381
rt_EI['Toluene'] = 2215.94723
rt_EI["Trichlorofluoromethane"]= 1597.372393
rt_EI["Dichlorodifluoromethane"] = 1468.148516
rt_EI["Tetrachloromethane"] = 1740.769826

dict_rt_EI_for_file['230311.0125.tank'] = rt_EI

#11
rt_EI['Benzene'] = 1908.479381
rt_EI['Toluene'] = 2215.94723
rt_EI["Trichlorofluoromethane"]= 1597.372393
rt_EI["Dichlorodifluoromethane"] = 1468.148516
rt_EI["Tetrachloromethane"] = 1740.769826

dict_rt_EI_for_file['230311.0555.tank'] = rt_EI

#14
rt_EI['Benzene'] = 1908.479381
rt_EI['Toluene'] = 2215.94723
rt_EI["Trichlorofluoromethane"]= 1597.372393
rt_EI["Dichlorodifluoromethane"] = 1468.148516
rt_EI["Tetrachloromethane"] = 1740.769826

dict_rt_EI_for_file['230315.0350.tank'] = rt_EI

# 
dict_rt_EI_for_file['230225.1804.std'] = rt_EI

directory_mass_cal = r"C:\Users\kaho\Desktop\data\data_Empa\mass_cal"
if not os.path.exists(directory_mass_cal):
    os.makedirs(directory_mass_cal)



for filename, rt_data in dict_rt_EI_for_file.items():
    with open(os.path.join(directory_mass_cal, f"{filename}.calibrants.EI.txt"), 'w') as file:
        for substance, retention_time in rt_data.items():
            file.write(f"{substance} = {retention_time}\n")