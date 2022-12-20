import numpy as np
import matplotlib.pyplot as plt

from itertools import product
from .io_functions import check_key,print_keys,copy_single_run

import matplotlib.pyplot as plt

from itertools import product
from .io_functions import check_key,print_keys,copy_single_run
from .vis_functions import *
from .fit_functions import gaussian,loggaussian,gaussian_train,loggaussian_train

def generate_cut_array(my_runs):
    """
    This function generates an array of bool = True with length = NEvts. If cuts are applied and then you run this function, it resets the cuts.
    """
    for run, ch in product(my_runs["NRun"], my_runs["NChannel"]):    
        for key in my_runs[run][ch].keys():
            if key.find("ADC") > 0:
                ADC_key = key
                # print(key)
        my_runs[run][ch]["MyCuts"] = np.ones(len(my_runs[run][ch][ADC_key]),dtype=bool)

def cut_min_max(my_runs, keys, limits, ranges = [0,0]):
    """
    This a fuction for cuts of min - max values. It takes a variable(s) and checks whether its value is between the specified limits.
    VARIABLES:
        - keys: a LIST of variables you want to constrain
        - limits: a DICTIONARY with same keys than variable "keys" and a list of the min and max values you want.
        - ranges: a LIST with the range where we want to check the key value. If [0,0] it uses the whole window. Time in sec.
    Important! Each key works independently. If one key gives True and the other False, it remains False.
    Example: keys = ["PeakAmp", "PeakTime"], limits = {"PeakAmp": [20,50], "PeakTime": [4e-6, 5e-6]}
    """
    for run, ch, key in product(my_runs["NRun"], my_runs["NChannel"], keys):
        if check_key(my_runs[run][ch], "MyCuts") == True:
            if check_key(my_runs[run][ch], key) == True:
                if ranges[0]==0 and ranges[1]==0:
                    for i in range(len(my_runs[run][ch][key])):
                        if limits[key][0] <= my_runs[run][ch][key][i] <= limits[key][1]:
                            continue    
                        else: my_runs[run][ch]["MyCuts"][i] = False
                else:
                    i_idx = int(np.round(ranges[0]/my_runs[run][ch]["Sampling"])); f_idx = int(np.round(ranges[1]/my_runs[run][ch]["Sampling"]))
                    for i in range(i_idx,f_idx+1):
                        if limits[key][0] <= my_runs[run][ch][key][i] <= limits[key][1]:
                            continue    
                        else: my_runs[run][ch]["MyCuts"][i] = False
            else: print(key," does not exist in my_runs!")
        else: print("Run generate_cut_array")

def cut_min_max_sim(my_runs, keys, limits):
    """
    This a fuction for cuts of min - max values. It takes a variable(s) and checks whether its value is between the specified limits.
    VARIABLES:
        - keys: a LIST of variables you want to constrain at the same time
        - limits: a DICTIONARY with same keys than variable "keys" and a list of the min and max values you want.
    Important! Keys are related, so all keys must be False to cut the event. If any of the conditions is True, the event is not cutted.
    Example: keys = ["PeakAmp"], limits = {"PeakAmp": [20,50]}
    """
    for run, ch in product(my_runs["NRun"], my_runs["NChannel"]):
        if check_key(my_runs[run][ch], "MyCuts") == True:
            for i in range(len(my_runs[run][ch][keys[0]])):
                for j in range(len(keys)):
                    if check_key(my_runs[run][ch], keys[j]) == True:
                        if limits[keys[j]][0] <= my_runs[run][ch][keys[j]][i] <= limits[keys[j]][1]:
                            my_runs[run][ch]["MyCuts"][i] = True
                            break
                        else: my_runs[run][ch]["MyCuts"][i] = False
                    else: print(keys," does not exist in my_runs!")
        else: print("Run generate_cut_array")

# def cut_std(my_runs, keys, limits):
#     for run, ch, key in product(my_runs["NRun"], my_runs["NChannel"], keys):
#         data = my_runs[run][ch][key]
#         ypbot = np.percentile(data, 0.1); yptop = np.percentile(data, 0.99)
#         ypad = 0.2*(yptop - ypbot)
#         ymin = ypbot - ypad; ymax = yptop + ypad
#         data = [i for i in data if ymin<i<ymax]

# def cut_std(my_runs, keys, limits):
#     for run, ch, key in product(my_runs["NRun"], my_runs["NChannel"], keys):
#         data = my_runs[run][ch][key]
#         ypbot = np.percentile(data, 0.1); yptop = np.percentile(data, 0.99)
#         ypad = 0.2*(yptop - ypbot)
#         ymin = ypbot - ypad; ymax = yptop + ypad
#         data = [i for i in data if ymin<i<ymax]
#         for run, ch, key in product(my_runs["NRun"], my_runs["NChannel"], keys):
#             if check_key(my_runs[run][ch], "MyCuts") == True:
#                 if check_key(my_runs[run][ch], key) == True:
#                     for i in range(len(my_runs[run][ch][key])):
#                         if limits[key][0] <= my_runs[run][ch][key][i] <= limits[key][1]:
#                             continue
#                         else: my_runs[run][ch]["MyCuts"][i] = False
#                 else: print(key," does not exist in my_runs!")
#             else: print("Run generate_cut_array")

def cut_lin_rel(my_runs, keys):
    for run, ch in product(my_runs["NRun"], my_runs["NChannel"]):
        if check_key(my_runs[run][ch], "MyCuts") == True:
                for j in range(len(keys)):
                    if check_key(my_runs[run][ch], keys[j]) == True:
                        continue
                    else: print("IAAA"); break
                plt.scatter(my_runs[run][ch][keys[0]], my_runs[run] [ch][keys[1]])
                # if limits[keys[j]][0] <= my_runs[run][ch][keys[j]][i] <= limits[keys[j]][1]:
                #     my_runs[run][ch]["MyCuts"][i] = True
                #     break
                # else: my_runs[run][ch]["MyCuts"][i] = False
        else: print("Run generate_cut_array")
