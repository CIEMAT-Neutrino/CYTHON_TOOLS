import sys
sys.path.insert(0, '../')

from lib.io_functions import load_npy,load_average_npy
from lib.dec_functions import deconvolve
from lib.dec_functions import deconvolve

run = 26
dec_run = 10
ch = 6

# my_runs = load_npy([run],[ch])
my_runs = load_average_npy([run],[ch]) # Activate in case deconvolution of average wvf wants to be performed
dec_runs = load_average_npy([dec_run],[ch])

OPT = {
    # "AVE":"ADC",
    "AVE":"AvWvf_threshold", # Activate in case deconvolution of average wvf wants to be performed
    "FIX_EXP":True,
    "LOGY":False,
    "FOCUS":False,
    "SHOW_F_SIGNAL":True,
    "SHOW_F_GAUSS":True,
    "SHOW_F_WIENER":True,
    "SHOW_F_DEC":True,
    "TRIMM": 0,
    "AUTO_TRIMM":True,
    "REVERSE":True,
    "SMOOTH":0
    }

# print(dec_runs[dec_run][ch].keys())
det_response = dec_runs[dec_run][ch]["AvWvf"]

deconvolve(my_runs,det_response,OPT)