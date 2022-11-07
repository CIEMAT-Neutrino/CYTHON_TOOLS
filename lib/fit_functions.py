import numpy as np
import scipy as sc
import matplotlib.pyplot as plt

from .io_functions import load_npy,load_analysis_npy
from .io_functions import check_key
from scipy.optimize import curve_fit
from scipy.special import erf

from itertools import product

def gaussian_train(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        amp = params[i+1]
        wid = params[i+2]
        y = y + amp * np.exp( -((x - ctr)/wid)**2)
    return y

def loggaussian_train(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        amp = params[i+1]
        wid = params[i+2]
        y = y + amp * np.exp( -((x - ctr)/wid)**2)
    return np.log10(y)

def gaussian(x, height, center, width):
    return height*np.exp(-(x - center)**2/(2*width**2))

def loggaussian(x, height, center, width):
    return np.log10(height*np.exp(-(x - center)**2/(2*width**2)))

def func(T,A,SIGMA,TAU,T0):
    return (2*A/TAU)*np.exp((SIGMA/(np.sqrt(2)*TAU))**2-(np.array(T)-T0)/TAU)*(1-erf((SIGMA**2-TAU*(np.array(T)-T0))/(np.sqrt(2)*SIGMA*TAU)))

def func2(T,P,A1,SIGMA,TAU1,T0,A2,TAU2):
    return P + func(T,A1,SIGMA,TAU1,T0) + func(T,A2,SIGMA,TAU2,T0)

def logfunc2(T,P,A1,SIGMA,TAU1,T0,A2,TAU2):
    return np.log(P + func(T,A1,SIGMA,TAU1,T0) + func(T,A2,SIGMA,TAU2,T0))

def logfunc3(T,P,A1,SIGMA,TAU1,T0,A2,TAU2,A3,TAU3):
    return np.log(P + func(T,A1,SIGMA,TAU1,T0) + func(T,A2,SIGMA,TAU2,T0) + func(T,A3,SIGMA,TAU3,T0))

def func3(T,P,A1,SIGMA,TAU1,T0,A2,TAU2,A3,TAU3):
    return P + func(T,A1,SIGMA,TAU1,T0) + func(T,A2,SIGMA,TAU2,T0) + func(T,A3,SIGMA,TAU3,T0)

def scfunc(T,A,B,C,D,E,F):
    return (A*np.exp(-(T-C)/B)/np.power(2*np.pi,0.5)*np.exp(-D**2/(B**2)))*(1-erf(((C-T)/D+D/B)/np.power(2,0.5)))+(E*np.exp(-(T-C)/F)/np.power(2*np.pi,0.5)*np.exp(-D**2/(F**2)))*(1-erf(((C-T)/D+D/F)/np.power(2,0.5)))

def sipm_fit(RAW,RAW_X,FIT_RANGE,OPT):
    MAX = np.argmax(RAW)
    thrld = 1e-1
    BUFFER1 = FIT_RANGE[0]
    BUFFER2 = FIT_RANGE[1]

    OPT["CUT_NEGATIVE"] = True
    popt1,perr1 = peak_fit(RAW, RAW_X,BUFFER1,OPT)

    p = np.mean(RAW[:MAX-BUFFER1])
    a1 = 2e-5; a1_low = 1e-8;  a1_high = 9e-2                    
    a2 = 2e-5; a2_low = 1e-8; a2_high = 9e-2                    
    a3 = 2e-5; a3_low = 1e-8; a3_high = 9e-2
    tau1 = 9e-8; tau1_low = 6e-9; tau1_high = 1e-7
    tau2 = 9e-7; tau2_low = tau1_high; tau2_high = 1e-6
    tau3 = 9e-6; tau3_low = tau2_high; tau3_high = 1e-5
    sigma2 = popt1[1]*10; sigma2_low = popt1[1]; sigma2_high = popt1[1]*100

    # USING VALUES FROM FIRST FIT PERFORM SECONG FIT FOR THE SLOW COMPONENT
    bounds2 = ([sigma2_low,a2_low,tau2_low,a3_low,tau3_low],[sigma2_high,a2_high,tau2_high,a3_high,tau3_high])
    initial2 = (sigma2,a2,tau2,a3,tau3)
    labels2 = ["SIGM","AMP2","TAU2","AMP3","TAU3"]
    popt, pcov = curve_fit(lambda T,SIGMA2,A2,TAU2,A3,TAU3: logfunc3(T,p,a1,SIGMA2,tau1,popt1[3],A2,TAU2,A3,TAU3),RAW_X[MAX-BUFFER1:MAX+BUFFER2],np.log(RAW[MAX-BUFFER1:MAX+BUFFER2]),p0 = initial2, bounds = bounds2, method = "trf")
    perr2 = np.sqrt(np.diag(pcov))

    sigma2 = popt[0];a2 = popt[1];tau2 = popt[2];a3 = popt[3];tau3 = popt[4]
    param = [p,a1,sigma2,tau1,popt1[3],a2,tau2,a3,tau3]
    if check_key(OPT,"SHOW") == True and OPT["SHOW"] == True:
        print("\n--- SECOND FIT VALUES (SLOW) ---")
        for i in range(len(initial2)):
            print("%s: \t%.2E \u00B1 %.2E"%(labels2[i],popt[i],perr2[i]))
        print("--------------------------------\n")

        # CHECK FIRST FIT
        plt.rcParams['figure.figsize'] = [16, 8]
        plt.subplot(1,2,1)
        plt.title("First fit to determine peak")
        plt.plot(RAW_X,RAW,label="RAW")
        plt.plot(RAW_X[MAX-BUFFER1:MAX+int(BUFFER1/2)],func(RAW_X[MAX-BUFFER1:MAX+int(BUFFER1/2)],*popt1),label="FIT")
        # plt.axvline(RAW_X[-buffer2],ls = "--",c = "k")
        plt.xlabel("Time in [s]"); plt.ylabel("ADC Counts")
        if check_key(OPT,"LOGY") != False: plt.semilogy();plt.ylim(thrld,RAW[MAX]*1.1)
        plt.legend()

        plt.subplot(1,2,2)
        plt.title("Second fit with full wvf")
        plt.plot(RAW_X,RAW,zorder=0,c="tab:blue",label="RAW")
        plt.plot(RAW_X[MAX-BUFFER1:MAX+BUFFER2],func3(RAW_X[MAX-BUFFER1:MAX+BUFFER2],*param),c="tab:orange",label="FIT")
        plt.plot(RAW_X,func3(RAW_X,*param),c="tab:green",label="FIT_FULL_LENGHT")
        plt.xlabel("Time in [s]"); plt.ylabel("ADC Counts")
        # plt.axvline(RAW_X[-buffer2],ls = "--",c = "k")
        if check_key(OPT,"LOGY") != False: plt.semilogy();plt.ylim(thrld,RAW[MAX]*1.1)
        plt.legend()
        while not plt.waitforbuttonpress(-1): pass
        plt.clf()

    aux = func3(RAW_X,*param)
    plt.ioff()
    return aux

def scint_fit(RAW,RAW_X,FIT_RANGE,OPT):        
    next_plot = False
    MAX = np.argmax(RAW)
    thrld = 1e-10
    BUFFER1 = FIT_RANGE[0]
    BUFFER2 = FIT_RANGE[1]

    OPT["CUT_NEGATIVE"] = True
    popt1,perr1 = peak_fit(RAW, RAW_X,BUFFER1,OPT)

    # USING VALUES FROM FIRST FIT PERFORM SECONG FIT FOR THE SLOW COMPONENT
    p = np.mean(RAW[:MAX-BUFFER1])
    p = 1e-1
    a1 = 5e-5; a1_low = 1e-6;  a1_high = 9e-3                               
    a3 = 5e-5; a3_low = 1e-7; a3_high = 9e-4
    tau1 = 1e-8; tau1_low = 6e-9; tau1_high = 9e-8
    tau3 = 5e-7; tau3_low = tau1_high; tau3_high = 9e-6
    sigma2 = popt1[1]; sigma2_low = popt1[1]*0.1; sigma2_high = popt1[1]*10
    bounds2 = ([sigma2_low,a3_low,tau3_low],[sigma2_high,a3_high,tau3_high])
    initial2 = (sigma2,a3,tau3)
    labels2 = ["SIGM","AMP2","TAU2"]
    popt2, pcov2 = curve_fit(lambda T,SIGMA2,A3,TAU3: logfunc2(T,p,popt1[0],SIGMA2,popt1[2],popt1[3],A3,TAU3),RAW_X[MAX-BUFFER1:MAX+BUFFER2],np.log(RAW[MAX-BUFFER1:MAX+BUFFER2]),p0 = initial2, bounds = bounds2, method = "trf")
    perr2 = np.sqrt(np.diag(pcov2))
    sigma2 = popt2[0];a3 = popt2[1];tau3 = popt2[2]
    param = [p,popt1[0],sigma2,popt1[2],popt1[3],a3,tau3]
    
    if check_key(OPT,"SHOW") == True and OPT["SHOW"] == True:
        print("\n--- SECOND FIT VALUES (SLOW) ---")
        for i in range(len(initial2)):
            print("%s: \t%.2E \u00B1 %.2E"%(labels2[i],popt2[i],perr2[i]))
        print("--------------------------------\n")

        print("SHOW key not included in OPT")
        # CHECK FIRST FIT
        plt.rcParams['figure.figsize'] = [16, 8]
        plt.subplot(1,2,1)
        plt.title("First fit to determine peak")
        plt.plot(RAW_X,RAW,label="RAW")
        plt.plot(RAW_X[MAX-BUFFER1:MAX+int(BUFFER1/2)],func(RAW_X[MAX-BUFFER1:MAX+int(BUFFER1/2)],*popt1),label="FIT")
        # plt.axvline(RAW_X[-buffer2],ls = "--",c = "k")
        plt.xlabel("Time in [s]"); plt.ylabel("ADC Counts")
        if check_key(OPT,"LOGY") == True and OPT["LOGY"] == True: plt.semilogy();plt.ylim(thrld,RAW[MAX]*1.1)
        plt.legend()
        plt.subplot(1,2,2)
        plt.title("Second fit with full wvf")
        plt.plot(RAW_X,RAW,zorder=0,c="tab:blue",label="RAW")
        # plt.plot(RAW_X[MAX-BUFFER:buffer2],func2sigma(RAW_X[MAX-BUFFER:buffer2],*param),c="tab:orange",label="FIT")
        plt.plot(RAW_X[MAX-BUFFER1:MAX+BUFFER2],func2(RAW_X[MAX-BUFFER1:MAX+BUFFER2],*param),c="tab:orange",label="FIT")
        # plt.plot(RAW_X,func2sigma(RAW_X,*param),c="tab:green",label="FIT_FULL_LENGHT")
        plt.xlabel("Time in [s]"); plt.ylabel("ADC Counts")
        plt.axvline(RAW_X[MAX+BUFFER2],ls = "--",c = "k")
        if check_key(OPT,"LOGY") == True and OPT["LOGY"] == True: plt.semilogy();plt.ylim(thrld,RAW[MAX]*1.1)
        plt.legend()
        while not plt.waitforbuttonpress(-1): pass
        plt.clf()
    aux = func2(RAW_X,*param)
    return aux

def sc_fit(RAW,RAW_X,FIT_RANGE,OPT):
    plt.ion()
    next_plot = False
    plt.rcParams['figure.figsize'] = [8, 8]
    FIT_RAW_X = np.arange(len(RAW))
    MAX = np.argmax(RAW)

    # USING VALUES FROM FIRST FIT PERFORM SECONG FIT FOR THE SLOW COMPONENT
    t0 = np.argmax(RAW)
    initial = (1500,150,t0,8,-700,300)
    labels = ["AMP","TAU1","T0","SIGMA","AMP2","TAU2"]

    try:
        popt, pcov = curve_fit(scfunc,FIT_RAW_X,RAW,p0 = initial, method = "trf")
        perr = np.sqrt(np.diag(pcov))
    except:
        print("Fit did not succeed")
        popt = initial
        perr = np.zeros(len(initial))

    if check_key(OPT, "SHOW") == True and OPT["SHOW"] == True:
        print("\n--- FIT VALUES (SLOW) ---")
        for i in range(len(initial)):
            print("%s: \t%.2E \u00B1 %.2E"%(labels[i],popt[i],perr[i]))
        print("--------------------------------\n")

        plt.title("Fit with full wvf")
        plt.plot(RAW_X,RAW,zorder=0,c="tab:blue",label="RAW")
        plt.plot(RAW_X,scfunc(FIT_RAW_X,*popt),"tab:orange",label="FIT")
        plt.xlabel("Time in [s]"); plt.ylabel("ADC Counts")
        if check_key(OPT,"LOGY") == True and OPT["LOGY"] == True: 
            plt.semilogy()
            plt.ylim(thrld,RAW[MAX]*1.1)
        plt.legend()
        while not plt.waitforbuttonpress(-1): pass
        plt.clf()
    aux = scfunc(FIT_RAW_X,*popt)
    return aux

def peak_fit(FIT_RAW,RAW_X,BUFFER,OPT):
    MAX = np.argmax(FIT_RAW)
    
    if check_key(OPT, "CUT_NEGATIVE") == True and OPT["CUT_NEGATIVE"] == True:
        for i in range(len(FIT_RAW)):
            if FIT_RAW[i] <= 1e-10:
                FIT_RAW[i] = 1e-10
            if np.isnan(FIT_RAW[i]):
                FIT_RAW[i] = 1e-10

    guess_t0 = RAW_X[np.argmax(FIT_RAW)-10]
    p = np.mean(FIT_RAW[:MAX-BUFFER])

    t0 = guess_t0; t0_low = guess_t0*0.02; t0_high = guess_t0*50
    sigma = 2e-8; sigma_low = 6e-9; sigma_high = 9e-8
    a1 = 2e-5; a1_low = 1e-8;  a1_high = 9e-2                    
    tau1 = 9e-8; tau1_low = 6e-9; tau1_high = 1e-7
    bounds = ([a1_low,sigma_low,tau1_low,t0_low],[a1_high,sigma_high,tau1_high,t0_high])
    initial = (a1,sigma,tau1,t0)
    labels = ["AMP1","SIG1","TAU1","TIME"]

    # FIT PEAK
    try:
        popt, pcov = curve_fit(func,RAW_X[MAX-BUFFER:MAX+int(BUFFER/2)],FIT_RAW[MAX-BUFFER:MAX+int(BUFFER/2)],p0 = initial, bounds = bounds, method = "trf")
        perr = np.sqrt(np.diag(pcov))
    except:
        print("Peak fit could not be performed")
        popt = initial
        perr = np.zeros(len(initial))

    # PRINT FIRST FIT VALUE
    if check_key(OPT,"SHOW") == True and OPT["SHOW"] == True:
        print("\n--- FISRT FIT VALUES (FAST) ---")
        for i in range(len(initial)):
            print("%s: \t%.2E \u00B1 %.2E"%(labels[i],popt[i],perr[i]))
        print("-------------------------------")

    # EXPORT FIT PARAMETERS
    a1 = popt[0];sigma = popt[1];tau1 = popt[2];t0 = popt[3]

    return popt,perr

def fit_wvfs(my_runs,signal_type,FIT_RANGE,OPT,PATH="../data/fit/"):
    
    try:
        ana_runs = load_analysis_npy(my_runs["N_runs"],my_runs["N_channels"])
    except:
        print("Events have not been processed")
    aux = dict()
    for run,ch in product(my_runs["N_runs"],my_runs["N_channels"]):
        if check_key(OPT, "AVE") == True and (OPT["AVE"] == "AvWvf" or OPT["AVE"] == "AvWvf_peak" or OPT["AVE"] == "AvWvf_threshold" or OPT["AVE"] == "Deconvolution"): 
            RAW = [my_runs[run][ch][OPT["AVE"]]]
            LOOP = 1
        else:
            RAW = ana_runs[run][ch]["P_channel"]*((my_runs[run][ch]["ADC"].T-ana_runs[run][ch]["Ped_mean"]).T)
            LOOP = len(my_runs[run][ch]["ADC"])
        
        RAW_X = 4e-9*np.arange(len(RAW[0]))
        for i in range(LOOP):
            print("Fitting wvf ",i)
            if signal_type == "SiPM":  fit = sipm_fit(RAW[i],RAW_X,FIT_RANGE,OPT)
            if signal_type == "SCINT": fit = scint_fit(RAW[i],RAW_X,FIT_RANGE,OPT)
            if signal_type == "SC":    fit = sc_fit(RAW[i],RAW_X,FIT_RANGE,OPT)
            aux[i] = fit
        
        plt.ioff()
        my_runs[run][ch][signal_type] = aux
        aux_path=PATH+"Fit_run"+str(run).zfill(2)+"_ch"+str(ch)+".npy"
        
        try:
            del my_runs[run][ch]["ADC"]
        except:
            print("'ADC' branch has already been deleted")

        np.save(aux_path,my_runs[run][ch])
        print("Saved data in:" , aux_path)