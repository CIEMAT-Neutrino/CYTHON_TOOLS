import sys
sys.path.insert(0, '../')

import matplotlib.pyplot as plt
import numpy as np

from .io_functions import load_npy,check_key
from itertools import product

def vis_npy(run,ch,KEY,OPT):
    
    norm_raw = 1
    norm_ave = 1
    buffer = 100
    
    RUN = load_npy([run],[ch])

    try:
        ANA_RUN = load_npy([run],[ch],"Analysis_","../data/ana/")
        ANA_BOOL = True
    except:
        ANA_BOOL = False
        print("Run has not been processed!")

    try:
        AVE_RUN = load_npy([run],[ch],"Average_","../data/ave/")
        ANA_BOOL = True
    except:
        AVE_BOOL = False
        print("Run has not been averaged!")
    
    sampling = RUN[run][ch]["Sampling"]
    sampling = 4e-9
    try:
        sampling = ANA_RUN[run][ch]["Sampling"]
    except:
        print("Showing sampling extracted from raw data")

    plt.ion()
    next_plot = False
    
    for i in range(len(RUN[run][ch]["ADC"])):
            
        plt.xlabel("Time in [s]")
        plt.ylabel("ADC Counts")
        plt.grid(True, alpha = 0.7)
        
        if (KEY == "ADC"):
            min = np.argmin(RUN[run][ch][KEY][i])
            RAW = RUN[run][ch][KEY][i]
            PED = np.mean(RUN[run][ch][KEY][i][:min-buffer])
            STD = np.std(RUN[run][ch][KEY][i][:min-buffer])
        
        elif(KEY == "Ana_ADC"):
            min = np.argmax(ANA_RUN[run][ch][KEY][i])
            RAW = ANA_RUN[run][ch][KEY][i]
            PED = 0    
            STD = ANA_RUN[run][ch]["Ped_STD"][i]

        if OPT["NORM"] == True and OPT["NORM"] == True:
            norm_raw = np.max(RAW)
            RAW = RAW/np.max(RAW)
        plt.plot(sampling*np.arange(len(RAW)),RAW,label="RAW_WVF")             

        if check_key(OPT, "SHOW_AVE") == True:   
            
            try:
                AVE_KEY = OPT["SHOW_AVE"]
                AVE = AVE_RUN[run][ch][AVE_KEY][0]
                if OPT["NORM"] == True and OPT["NORM"] == True:
                    AVE = AVE/np.max(AVE)
                plt.plot(AVE_RUN[run][ch]["Sampling"]*np.arange(len(AVE)),AVE,alpha=.5,label="AVE_WVF_%s"%AVE_KEY)             

            except:
                print("Run has not been averaged!")
                      
        
        if OPT["LOGY"] == True:
            plt.semilogy()

        plt.plot(4e-9*np.array([min-buffer,min-buffer]),np.array([PED+5*STD,PED-5*STD])/norm_raw,c="red",lw=2.)
        plt.axhline((PED)/norm_raw,c="k",alpha=.5)
        plt.axhline((PED+STD)/norm_raw,c="k",alpha=.5,ls="--")
        plt.axhline((PED-STD)/norm_raw,c="k",alpha=.5,ls="--")
        plt.legend()

        if OPT["SHOW_PARAM"] == True:
            print("Event Number {} from RUN_{} CH_{} ({})".format(i,run,ch,ANA_RUN[run][ch]["Label"]))
            print("Sampling: {:.0E}".format(ANA_RUN[run][ch]["Sampling"]))
            print("Pedestal mean: {:.2E}".format(ANA_RUN[run][ch]["Ped_mean"][i]))
            print("Pedestal STD: {:.4f}".format(ANA_RUN[run][ch]["Ped_STD"][i]))
            print("Pedestal min: {:.4f}\t Pedestal max {:.4f}".format(ANA_RUN[run][ch]["Ped_min"][i],ANA_RUN[run][ch]["Ped_max"][i]))
            print("Pedestal time limit: {:.4E}".format(4e-9*(min-buffer)))
            print("Max Peak Amplitude: {:.4f}".format(ANA_RUN[run][ch]["Peak_amp"][i]))
            print("Max Peak Time: {:.2E}".format(ANA_RUN[run][ch]["Peak_time"][i]*ANA_RUN[run][ch]["Sampling"]))
            print("Charge: {:.2E}\n".format(ANA_RUN[run][ch]["AVE_INT_LIMITS"][i]))

        while not plt.waitforbuttonpress(-1): pass

        plt.clf()

    plt.ioff()

def vis_ave_npy(run,ch,CUTS,OPT,PATH = "../data/raw"):
    aux_ADC=ana_runs[run][ch]["P_channel"]*((my_runs[run][ch]["ADC"].T-ana_runs[run][ch]["Ped_mean"]).T)
    my_runs[run][ch]["AvWvf"] = np.mean(aux_ADC,axis=0)
    return my_runs[run][ch]["AvWvf"]

def vis_var_hist(my_run,KEY,w=1e-4):
    plt.ion()
    for run, ch, key in product(my_run["N_runs"],my_run["N_channels"],KEY):
        w=abs(np.max(my_run[run][ch][key])-np.min(my_run[run][ch][key]))*w
        data = []
        if key == "Peak_amp":
            data = my_run[run][ch][key]
            max_amp = np.max(data)
            binning = int(max_amp)+1
        elif key == "Peak_time":
            data = 4e-9*my_run[run][ch][key]
            binning = int(my_run[run][ch]["NBins_wvf"]+1)
        else:
            data = my_run[run][ch][key]
            data=sorted(data)
            out_threshold= 2.0*np.std(data+[-a for a in data]) # Repeat distribution in negative axis to compute std
            data=[i for i in data if -out_threshold<i<out_threshold]
            binning = np.arange(min(data), max(data) + w, w)
            if len(binning) < 100: binning = 100

        plt.hist(data,binning)
        while not plt.waitforbuttonpress(-1): pass
        plt.clf()
    plt.ioff()