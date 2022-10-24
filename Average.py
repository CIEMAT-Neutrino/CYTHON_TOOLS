from lib.my_functions import load_npy
from lib.wvf_functions import average_wvf

N_runs     =[10,22,26]     
N_channels =[0,1,4,6]       
Pl         =[-1,-1,-1,-1]   #polarity
P_channels ={}
for ch,pl in zip(N_channels,Pl): P_channels[ch]=pl

L_channels  =["SiPM1","SiPM2","PMT","SuperCell"]
RUNS = load_npy(N_runs, N_channels,P_channels,"data/")

average_wvf(RUNS)