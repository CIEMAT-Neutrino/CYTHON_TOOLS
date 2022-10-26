# import copy
import numpy as np
from .io_functions import load_analysis_npy, load_npy
from itertools import product

def average_wvfs(my_runs,PATH="../data/ave/",threshold=50):
    
    try:
        ana_runs = load_analysis_npy(my_runs["N_runs"],my_runs["N_channels"])
    except:
        print("EVENTS HAVE NOT BEEN PROCESSED! Please run Process.py")
        return 0

    for run,ch in product(my_runs["N_runs"],my_runs["N_channels"]):
        
        # No centering
        aux_ADC=ana_runs[run][ch]["P_channel"]*((my_runs[run][ch]["ADC"].T-ana_runs[run][ch]["Ped_mean"]).T);

        my_runs[run][ch]["AvWvf"] = np.mean(aux_ADC,axis=0)
        aux_path=PATH+"Average_run"+str(run).zfill(2)+"_ch"+str(ch)+".npy"
        
        # centering
        bin_ref         =np.argmax(aux_ADC[0]) #using the first peak as reference
        # bin_ref         =int(len(aux_ADC[0])/15) #10% of time window
        av_wvf_peak      =np.zeros(aux_ADC[0].shape)
        av_wvf_threshold =np.zeros(aux_ADC[0].shape)

        n_wvs   =len(aux_ADC)
        n_bins  =len(aux_ADC[0])
        for wvf in aux_ADC:
            bin_peak      = np.argmax(wvf)
            
            try:
                bin_threshold = np.argwhere(wvf>threshold)[0][0]
            except:
                #no good value found;
                bin_threshold=n_bins

            # Peak centering
            if bin_ref<bin_peak:
                av_wvf_peak[:bin_ref]                    += (wvf[(bin_peak-bin_ref):bin_peak]/n_wvs);
                av_wvf_peak[bin_ref:-(bin_peak-bin_ref)] += (wvf[bin_peak:]/n_wvs);
            else:
                av_wvf_peak[(bin_ref-bin_peak):bin_ref] += (wvf[:bin_peak]/n_wvs);
                av_wvf_peak[bin_ref:-1]                   += (wvf[bin_peak:-(bin_ref-bin_peak+1)]/n_wvs);
            
            # threshold centering
            if bin_ref<bin_threshold:
                av_wvf_threshold[:bin_ref]                    += (wvf[(bin_threshold-bin_ref):bin_threshold]/n_wvs);
                av_wvf_threshold[bin_ref:-(bin_threshold-bin_ref)] += (wvf[bin_threshold:]/n_wvs);
            else:
                av_wvf_threshold[(bin_ref-bin_threshold):bin_ref] += (wvf[:bin_threshold]/n_wvs);
                av_wvf_threshold[bin_ref:-1]                   += (wvf[bin_threshold:-(bin_ref-bin_threshold+1)]/n_wvs);
        
        
        my_runs[run][ch]["AvWvf_peak"]=av_wvf_peak;
        my_runs[run][ch]["AvWvf_threshold"]=av_wvf_threshold;
        del my_runs[run][ch]["ADC"];

        np.save(aux_path,my_runs[run][ch])
        print("Saved data in:" , aux_path)

def integrate(my_runs,PATH="../data/ana/"):
    
    try:
        ana_runs = load_analysis_npy(my_runs["N_runs"],my_runs["N_channels"])
    except:
        print("EVENTS HAVE NOT BEEN PROCESSED! Please run Process.py")
        return 0

    for run,ch in product(my_runs["N_runs"],my_runs["N_channels"]):
        for i in range(len(my_runs[run][ch]["ADC"])):
            RAW = ana_runs[run][ch]["P_channel"]*(my_runs[run][ch]["ADC"][i]-ana_runs[run][ch]["Ped_mean"][i])
            MAX = np.argmax(RAW)
            INT_I = 0
            INT_F = 0
            
            for j in range(len(RAW[MAX:])):
                if RAW[MAX+j] < 0:
                    INT_F = MAX+j
                    break
            for j in range(len(RAW[:MAX])):
                if RAW[MAX-j] < 0:
                    INT_I = MAX-j+1
                    break
            print(RAW[INT_I:INT_F])
            my_runs[run][ch]["Int"][i] = np.trapz(RAW[INT_I:INT_F],x=4e-9*np.arange(len(RAW[INT_I:INT_F])))
    # TODO: save to .npy file