import sys
sys.path.insert(0, '../')

from lib.io_functions import load_npy,load_average_npy
from lib.dec_functions import deconvolve
from lib.dec_functions import deconvolve

run = 26
dec_run = 10
ch = 6

my_runs = load_npy([run],[ch])
dec_runs = load_average_npy([dec_run],[ch])

OPT = {"LOGY":False, "FOCUS":True}
OPT = {"LOGY":False}

# print(dec_runs[10][0].keys())
# KERNEL = dec_runs["ADC"]

deconvolve(my_runs,dec_runs[dec_run][ch]["AvWvf_peak"],120,0,OPT)