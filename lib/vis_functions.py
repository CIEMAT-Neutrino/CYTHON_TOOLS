import sys
sys.path.insert(0, '../')

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.cm import viridis
import numpy as np
import keyboard
import math

from .io_functions import load_npy,check_key
from .cut_functions import *
from .cut_functions import *
from itertools import product

import scipy
from scipy.signal import find_peaks

from .fig_config import (
    add_grid,
    figure_features,
)  # <--- import customized functions

def vis_npy(my_run, keys, OPT = {}, evt_sel = -1, same_plot = False):
    """
    This function is a event visualizer. It plots individual events of a run, indicating the pedestal level, pedestal std and the pedestal calc limit.
    We can interact with the plot and pass through the events freely (go back, jump to a specific event...)
    VARIABLES:
        - my_run: run(s) we want to check
        - KEYS: choose between ADC or AnaADC to see raw (as get from ADC) or Analyzed events (starting in 0 counts), respectively. Type: List
        - OPT: several options that can be True or False. Type: List
            a) MICRO_SEC: if True we multiply Sampling by 1e6
            b) NORM: True if we want normalized waveforms
            c) LOGY: True if we want logarithmic y-axis
            d) SHOW_AVE: if computed and True, it will show average
            e) SHOW_PARAM: True if we want to check calculated parameters (pedestal, amplitude, charge...)
            f) CHARGE_KEY: if computed and True, it will show the parametre value
            g) PEAK_FINDER: True if we want to check how many peaks are
        - evt_sel: choose the events we want to see. If -1 all events are displayed, if 0 only uncutted events are displayed, if 1 only cutted events are displayed
        - same_plot: True if we want to plot different channels in the SAME plot
    """

    charge_key = "ChargeAveRange"
    if check_key(OPT, "CHARGE_KEY"): charge_key = OPT["CHARGE_KEY"]
    norm_ave = 1
    buffer = 100
    figure_features()

    plt.ion()
    ch_list = my_run["NChannel"]
    nch = len(my_run["NChannel"])
    axs = []
    if same_plot == False:
        if nch < 4:
            fig, ax = plt.subplots(nch ,1, figsize = (10,8))
            if nch == 1: axs.append(ax)
            else: axs = ax
        else:
            fig, ax = plt.subplots(2, math.ceil(nch/2), figsize = (10,8))
            axs = ax.T.flatten()
    else:
        fig, ax = plt.subplots(1 ,1, figsize = (8,6))
        axs = ax

    for run, key in product(my_run["NRun"],keys):
        idx = 0
        for i in range(len(my_run[run][ch_list[0]][key])):
            try:
                skip = 0
                for ch in ch_list:
                    if evt_sel == 0 and my_run[run][ch]["MyCuts"][idx] == False: skip = 1; break # To Skip Cutted events!!
                    if evt_sel == 1 and my_run[run][ch]["MyCuts"][idx] == True: skip = 1; break # To Get only Cutted events!!
                if skip == 1: idx = idx +1; continue
            except: pass
            
            fig.supxlabel(r'Time [s]')
            fig.supylabel("ADC Counts")
            min = []
            raw = []
            norm_raw = [1]*nch # Generates a list with the norm correction for std bar
            for j in range(nch):
                if (key == "ADC"):
                    min.append(np.argmin(my_run[run][ch_list[j]][key][idx]))
                    raw.append(my_run[run][ch_list[j]][key][idx])
                    ped = np.mean(my_run[run][ch_list[j]][key][idx][:min[j]-buffer])
                    std = np.std(my_run[run][ch_list[j]][key][idx][:min[j]-buffer])
                
                if(key == "AnaADC"):
                    min.append(np.argmax(my_run[run][ch_list[j]][key][idx]))
                    raw.append(my_run[run][ch_list[j]][key][idx])
                    ped = 0
                    std = my_run[run][ch_list[j]]["PedSTD"][idx]

                if OPT["NORM"] == True and OPT["NORM"] == True:
                    norm_raw[j] = (np.max(raw[j]))
                    raw[j] = raw[j]/np.max(raw[j])

                sampling = my_run[run][ch_list[j]]["Sampling"] # To reset the sampling to its initial value (could be improved)
                if check_key(OPT, "MICRO_SEC") == True and OPT["MICRO_SEC"]==True:
                    fig.supxlabel(r'Time [$\mu$s]')
                    my_run[run][ch_list[j]]["Sampling"] = my_run[run][ch_list[j]]["Sampling"]*1e6

                if same_plot == False:
                    if OPT["LOGY"] == True:
                        axs[j].semilogy()
                        std = 0 # It is ugly if we see this line in log plots
                    fig.tight_layout(h_pad=2) # We avoid small vertical space between plots            
                    axs[j].plot(my_run[run][ch_list[j]]["Sampling"]*np.arange(len(raw[j])),raw[j],label="RAW_WVF", drawstyle = "steps", alpha = 0.95, linewidth=1.2)
                    axs[j].grid(True, alpha = 0.7)
                    axs[j].plot(my_run[run][ch_list[j]]["Sampling"]*np.array([my_run[run][ch_list[j]]["PedLim"],my_run[run][ch_list[j]]["PedLim"]]),np.array([ped+4*std,ped-4*std])/norm_raw[j],c="red",lw=2., alpha = 0.8)
                    axs[j].axhline((ped)/norm_raw[j],c="k",alpha=.55)
                    axs[j].axhline((ped+std)/norm_raw[j],c="k",alpha=.5,ls="--"); axs[j].axhline((ped-std)/norm_raw[j],c="k",alpha=.5,ls="--")
                    axs[j].set_title("Run {} - Ch {} - Event Number {}".format(run,ch_list[j],idx),size = 14)
                    axs[j].xaxis.offsetText.set_fontsize(14) # Smaller fontsize for scientific notation
                    
                    if check_key(OPT, "SHOW_AVE") == True:   
                        try:
                            ave_key = OPT["SHOW_AVE"]
                            ave = my_run[run][ch_list[j]][ave_key][0]
                            if OPT["NORM"] == True and OPT["NORM"] == True:
                                ave = ave/np.max(ave)
                            axs[j].plot(my_run[run][ch_list[j]]["Sampling"]*np.arange(len(ave)),ave,alpha=.5,label="AVE_WVF_%s"%ave_key)             
                        except:
                            print("Run has not been averaged!")

                    if check_key(OPT, "LEGEND") == True and OPT["LEGEND"]:
                        axs[j].legend()

                    if OPT["PEAK_FINDER"]:
                        # These parameters must be modified according to the run...
                        thresh = my_run[run][ch_list[j]]["PedMax"][idx]
                        wdth = 4
                        prom = 0.01
                        dist  = 40
                        axs[j].axhline(thresh,c="k", alpha=.6, ls = "dotted")
                        peak_idx, _ = find_peaks(raw[j], height = thresh, width = wdth, prominence = prom, distance=dist)
                        for p in peak_idx:
                            axs[j].scatter(my_run[run][ch_list[j]]["Sampling"]*p,raw[j][p],c="tab:red", alpha = 0.9)

                    try:
                        if my_run[run][ch_list[j]]["MyCuts"][idx] == False:
                            figure_features(tex = False)
                            axs[j].text(0.5,0.5, s = 'CUT', fontsize = 100, horizontalalignment='center',verticalalignment='center',
                                        transform = axs[j].transAxes, color = 'red', fontweight = "bold", alpha = 0.5)
                            figure_features()
                    except: pass
                
                else:
                    if OPT["LOGY"] == True:
                        axs.semilogy()
                        std = 0 # It is ugly if we see this line in log plots
                    axs.plot(my_run[run][ch_list[j]]["Sampling"]*np.arange(len(raw[j])),raw[j], drawstyle = "steps", alpha = 0.95, linewidth=1.2, label = "Ch {} ({})".format(ch_list[j],my_run[run][ch_list[j]]["Label"]))
                    axs.grid(True, alpha = 0.7)
                    axs.plot(my_run[run][ch_list[j]]["Sampling"]*np.array([my_run[run][ch_list[j]]["PedLim"],my_run[run][ch_list[j]]["PedLim"]]),np.array([ped+4*std,ped-4*std])/norm_raw[j],c="red",lw=2., alpha = 0.8)
                    # axs.axhline((ped)/norm_raw[j],c="k",alpha=.55)
                    # axs.axhline((ped+std)/norm_raw[j],c="k",alpha=.5,ls="--"); axs.axhline((ped-std)/norm_raw[j],c="k",alpha=.5,ls="--")
                    axs.set_title("Run {} - Event Number {}".format(run,idx),size = 14)
                    axs.xaxis.offsetText.set_fontsize(14)
                    
                    if check_key(OPT, "SHOW_AVE") == True:   
                        try:
                            ave_key = OPT["SHOW_AVE"]
                            ave = my_run[run][ch_list[j]][ave_key][0]
                            if OPT["NORM"] == True and OPT["NORM"] == True:
                                ave = ave/np.max(ave)
                            axs.plot(my_run[run][ch_list[j]]["Sampling"]*np.arange(len(ave)),ave,alpha=.5,label="AVE_WVF_%s"%ave_key)             
                        except:
                            print("Run has not been averaged!")

                    if check_key(OPT, "LEGEND") == True and OPT["LEGEND"]:
                        axs.legend()

                    if OPT["PEAK_FINDER"]:
                        # These parameters must be modified according to the run...
                        thresh = my_run[run][ch_list[j]]["PedMax"][idx]
                        wdth = 4
                        prom = 0.01
                        dist  = 40
                        axs.axhline(thresh,c="k", alpha=.6, ls = "dotted")
                        peak_idx, _ = find_peaks(raw[j], height = thresh, width = wdth, prominence = prom, distance=dist)
                        for p in peak_idx:
                            axs.scatter(my_run[run][ch_list[j]]["Sampling"]*p,raw[j][p],c="tab:red", alpha = 0.9)

                    try:
                        if my_run[run][ch_list[j]]["MyCuts"][idx] == False:
                            figure_features(tex = False)
                            axs.text(0.5,0.5, s = 'CUT', fontsize = 100, horizontalalignment='center',verticalalignment='center',
                                        transform = axs.transAxes, color = 'red', fontweight = "bold", alpha = 0.5)
                            figure_features()
                    except: pass
                    
                if OPT["SHOW_PARAM"] == True:
                    print('\033[1m' + "\nEvent Number {} from RUN_{} CH_{} ({})".format(idx,run,ch_list[j],my_run[run][ch_list[j]]["Label"]) + '\033[0m')
                    print("- Sampling: {:.0E}".format(sampling))
                    print("- Pedestal mean: {:.2E}".format(my_run[run][ch_list[j]]["PedMean"][idx]))
                    print("- Pedestal std: {:.4f}".format(my_run[run][ch_list[j]]["PedSTD"][idx]))
                    print("- Pedestal min: {:.4f}\t Pedestal max {:.4f}".format(my_run[run][ch_list[j]]["PedMin"][idx],my_run[run][ch_list[j]]["PedMax"][idx]))
                    print("- Pedestal time limit: {:.4E}".format(my_run[run][ch_list[j]]["Sampling"]*my_run[run][ch_list[j]]["PedLim"]))
                    print("- Max Peak Amplitude: {:.4f}".format(my_run[run][ch_list[j]]["PeakAmp"][idx]))
                    print("- Max Peak Time: {:.2E}".format(my_run[run][ch_list[j]]["PeakTime"][idx]*my_run[run][ch_list[j]]["Sampling"]))
                    try:
                        print("- Charge: {:.2E}".format(my_run[run][ch_list[j]][OPT["CHARGE_KEY"]][idx]))
                    except:
                        if check_key(OPT,"CHARGE_KEY"): print("- Charge: has not been computed for key %s!"%OPT["CHARGE_KEY"])
                        else: print("- Charge: default charge key has not been computed")
                    try:
                        print("- Peak_idx:",peak_idx*my_run[run][ch_list[j]]["Sampling"])
                    except:
                        if not check_key(OPT,"PEAK_FINDER"): print("")
                my_run[run][ch_list[j]]["Sampling"] = sampling    

            tecla = input("\nPress q to quit, p to save plot, r to go back, n to choose event or any key to continue: ")
            if tecla == "q":
                break
            elif tecla == "r":
                idx = idx-1
            elif tecla == "n":
                ev_num = int(input("Enter event number: "))
                idx = ev_num
                if idx > len(my_run[run][ch_list[j]][key]): idx = len(my_run[run][ch_list[j]][key])-1; print('\033[1m' + "\nBe careful! There are ", idx, "in total"); print('\033[0m')
            elif tecla == "p":
                fig.savefig('run{}_evt{}.png'.format(run,idx), dpi = 500)
                idx = idx+1
            else:
                idx = idx + 1
            if idx == len(my_run[run][ch_list[j]][key]): break
            try: [axs[j].clear() for j in range (nch)]
            except: axs.clear()
        try: [axs[j].clear() for j in range (nch)];
        except: axs.clear()
    plt.ioff()
    plt.clf()

def vis_var_hist(my_run, keys, percentile = [0.1, 99.9], OPT = {}):
    """
    This function takes the specified variables and makes histograms. The binning is fix to 600, so maybe it is not the appropriate.
    Outliers are taken into account with the percentile. It discards values below and above the indicated percetiles.
    It returns values of counts, bins and bars from the histogram to be used in other function.
    VARIABLES:
        - my_run: run(s) we want to check
        - keys: variables we want to plot as histograms. Type: List
            a) PeakAmp: histogram of max amplitudes of all events. The binning is 1 ADC. There are not outliers.
            b) PeakTime: histogram of times of the max amplitude in events. The binning is the double of the sampling. There are not outliers.
            c) Other variable: any other variable. Here we reject outliers.
        - percentile: percentile used for outliers removal
    WARNING! Maybe the binning stuff should be studied in more detail.
    """
    
    # keys is the variable that we want to plot
    w=1e-4 # w is related to the bin width in some way
    figure_features()

    plt.ion()
    all_counts = []
    all_bins = []
    all_bars = []
    for run, ch, key in product(my_run["NRun"],my_run["NChannel"],keys):
        # w = abs(np.max(my_run[run][ch][key])-np.min(my_run[run][ch][key]))*w
        data = []
        if key == "PeakAmp":
            data = my_run[run][ch][key]
            max_amp = np.max(data)
            binning = int(max_amp)+1
        elif key == "PeakTime":
            data = my_run[run][ch]["Sampling"]*my_run[run][ch][key]
            binning = int(my_run[run][ch]["NBinsWvf"]/2)
        else:
            data = my_run[run][ch][key]
            ypbot = np.percentile(data, percentile[0]); yptop = np.percentile(data, percentile[1])
            ypad = 0.2*(yptop - ypbot)
            ymin = ypbot - ypad; ymax = yptop + ypad
            data = [i for i in data if ymin<i<ymax]
            w = abs(np.max(data)-np.min(data))*w
            binning = 600 # FIXED VALUE UNTIL BETTER SOLUTION
            # binning = np.arange(min(data), max(data) + w, w)
            # binning = int((np.max(data)-np.min(data)) * len(data)**(1/3) / (3.49*np.std(data)))*5 # Visto por internete
            # print(w)
            # print(len(binning))
            # if len(binning) < 100: binning = 100
            # if len(binning) > 2000: binning = 1000
            # AALTERNATIVE WAY OF FILTERING OUTLIERS
            # data=sorted(data)
            # out_threshold= 2.0*np.std(data+[-a for a in data]) # Repeat distribution in negative axis to compute std
            # data=[i for i in data if -out_threshold<i<out_threshold]
            # binning = np.arange(min(data), max(data) + w, w)


        counts, bins, bars = plt.hist(data,binning) # , zorder = 2 f
        all_counts.append(counts)
        all_bins.append(bins)
        all_bars.append(bars)
        plt.grid(True, alpha = 0.7) # , zorder = 0 for grid behind hist
        plt.title("Run_{} Ch_{} - {} histogram".format(run,ch,key),size = 14)
        plt.xticks(size = 11); plt.yticks(size = 11)
        plt.xlabel(key+" (Units)", size = 11); plt.ylabel("Counts", size = 11)
        while not plt.waitforbuttonpress(-1): pass
        plt.clf()
    plt.ioff()
    return all_counts, all_bins, all_bars

def vis_persistence(my_run, OPT = {}):
    """
    This function plot the PERSISTENCE histogram of the given runs&ch.
    It perfoms a cut in 20<"PeakTime"(bins)<50 so that all the events not satisfying the condition are removed. 
    Binning is fixed (x=5000, y=1000) [study upgrade].
    X_data (time) and Y_data (waveforms) are deleted after the plot to save space.
    WARNING! flattening long arrays leads to MEMORY problems :/
    """

    plt.ion()
    for run, ch in product(my_run["NRun"],my_run["NChannel"]):

        generate_cut_array(my_run)
        cut_min_max(my_run, ["PeakTime"], {"PeakTime":[my_run[run][ch]["PedLim"]-20,my_run[run][ch]["PedLim"]+50]})

        data_flatten = my_run[run][ch]["AnaADC"][np.where(my_run[run][ch]["MyCuts"] == True)].flatten() #####
        time = my_run[run][ch]["Sampling"]*np.arange(len(my_run[run][ch]["AnaADC"][0]))
        time_flatten = np.array([time] * int(len(data_flatten)/len(time))).flatten()

        plt.hist2d(time_flatten,data_flatten,density=True,bins=[5000,1000], cmap = viridis, norm=LogNorm())

        plt.colorbar()
        plt.grid(True, alpha = 0.7) # , zorder = 0 for grid behind hist
        plt.title("Run_{} Ch_{} - Persistence".format(run,ch),size = 14)
        plt.xticks(size = 11); plt.yticks(size = 11)
        plt.xlabel("Time (s)", size = 11); plt.ylabel("Counts", size = 11)
        del data_flatten, time, time_flatten
        while not plt.waitforbuttonpress(-1): pass
        plt.clf()
    plt.ioff()
    plt.clf()

