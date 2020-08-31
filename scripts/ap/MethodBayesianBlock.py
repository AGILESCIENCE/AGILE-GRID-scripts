import numpy as np
from astropy.stats import bayesian_blocks
import matplotlib.pyplot as plt
import copy
from scipy import stats
import pandas as pd


binned_e1dq1 = np.fromfile('datasets/LC_3C454.3/E1dq1_21600s_emin100_emax10000_r4.ap.ap3', sep = ' ')
binned_e1dq1 = binned_e1dq1.reshape(len(binned_e1dq1)//27,27)
binned_e1dq1 = binned_e1dq1[2:,]

#### Fixed Exposre ####

binned_sim01_fixed = np.fromfile('datasets/LC_3C454.3/fixed_exp/E1dq1_21600s_emin100_emax10000_r4_run1.ap.sim.ap.ap3', sep = ' ')
binned_sim01_fixed = binned_sim01_fixed.reshape(len(binned_sim01_fixed)//27,27)
binned_sim01_fixed = binned_sim01_fixed[2:,]

binned_sim02_fixed = np.fromfile('datasets/LC_3C454.3/fixed_exp/E1dq1_21600s_emin100_emax10000_r4_run2.ap.sim.ap.ap3', sep = ' ')
binned_sim02_fixed = binned_sim02_fixed.reshape(len(binned_sim02_fixed)//27,27)
binned_sim02_fixed = binned_sim02_fixed[2:,]

binned_sim03_fixed = np.fromfile('datasets/LC_3C454.3/fixed_exp/E1dq1_21600s_emin100_emax10000_r4_run3.ap.sim.ap.ap3', sep = ' ')
binned_sim03_fixed = binned_sim03_fixed.reshape(len(binned_sim03_fixed)//27,27)
binned_sim03_fixed = binned_sim03_fixed[2:,]

binned_sim04_fixed = np.fromfile('datasets/LC_3C454.3/fixed_exp/E1dq1_21600s_emin100_emax10000_r4_run4.ap.sim.ap.ap3', sep = ' ')
binned_sim04_fixed = binned_sim04_fixed.reshape(len(binned_sim04_fixed)//27,27)
binned_sim04_fixed = binned_sim04_fixed[2:,]

binned_sim05_fixed = np.fromfile('datasets/LC_3C454.3/fixed_exp/E1dq1_21600s_emin100_emax10000_r4_run5.ap.sim.ap.ap3', sep = ' ')
binned_sim05_fixed = binned_sim05_fixed.reshape(len(binned_sim05_fixed)//27,27)
binned_sim05_fixed = binned_sim05_fixed[2:,]

binned_sim06_fixed = np.fromfile('datasets/LC_3C454.3/fixed_exp/E1dq1_21600s_emin100_emax10000_r4_run6.ap.sim.ap.ap3', sep = ' ')
binned_sim06_fixed = binned_sim06_fixed.reshape(len(binned_sim06_fixed)//27,27)
binned_sim06_fixed = binned_sim06_fixed[2:,]

binned_sim07_fixed = np.fromfile('datasets/LC_3C454.3/fixed_exp/E1dq1_21600s_emin100_emax10000_r4_run7.ap.sim.ap.ap3', sep = ' ')
binned_sim07_fixed = binned_sim07_fixed.reshape(len(binned_sim07_fixed)//27,27)
binned_sim07_fixed = binned_sim07_fixed[2:,]

binned_sim08_fixed = np.fromfile('datasets/LC_3C454.3/fixed_exp/E1dq1_21600s_emin100_emax10000_r4_run8.ap.sim.ap.ap3', sep = ' ')
binned_sim08_fixed = binned_sim08_fixed.reshape(len(binned_sim08_fixed)//27,27)
binned_sim08_fixed = binned_sim08_fixed[2:,]

binned_sim09_fixed = np.fromfile('datasets/LC_3C454.3/fixed_exp/E1dq1_21600s_emin100_emax10000_r4_run9.ap.sim.ap.ap3', sep = ' ')
binned_sim09_fixed = binned_sim09_fixed.reshape(len(binned_sim09_fixed)//27,27)
binned_sim09_fixed = binned_sim09_fixed[2:,]

binned_sim10_fixed = np.fromfile('datasets/LC_3C454.3/fixed_exp/E1dq1_21600s_emin100_emax10000_r4_run10.ap.sim.ap.ap3', sep = ' ')
binned_sim10_fixed = binned_sim10_fixed.reshape(len(binned_sim10_fixed)//27,27)
binned_sim10_fixed = binned_sim10_fixed[2:,]


# 0:tstart 
# 1:tstop 
# 2:exp[cm2 s]
# 3:cts  N_s + N_b 
# 4:normAB11 
# 5:normAB12 
# 6:normAB13 
# 7:normAB14 
# 8:normAB21 
# 9:normAB22 
#10:normAB23 
#11:normAB24 
#12:normAB11aa 
#13:normAB21aa 
#14:ratediffR1 
#15:ratediffR2 
#16:ratediffR3 
#17:ratediffR4  (N_s + N_b)/exp - N_b/exp
#18:ratediffR1AA 
#19:rate (compreso background) (N_s + N_b)/exp
#20:rate_error 
#21:flux_ratediffR4  
#22:flux_ratediffR4_error 
#23:TS # Test Statistico
#24:flux_rate 
#25:flux_rate_error 
#26:cts_expBKG4

# Significanza dei picchi
# real time analysis (un valore alla volta) OK


#### DATASET FIXED EXPOSURE ####
a0 = (binned_sim01_fixed[:,0],binned_sim02_fixed[:,0],binned_sim03_fixed[:,0],
      binned_sim04_fixed[:,0],binned_sim05_fixed[:,0],binned_sim06_fixed[:,0],
      binned_sim07_fixed[:,0],binned_sim08_fixed[:,0],binned_sim09_fixed[:,0],
      binned_sim10_fixed[:,0])
a19 = (binned_sim01_fixed[:,19],binned_sim02_fixed[:,19],binned_sim03_fixed[:,19],
      binned_sim04_fixed[:,19],binned_sim05_fixed[:,19],binned_sim06_fixed[:,19],
      binned_sim07_fixed[:,19],binned_sim08_fixed[:,19],binned_sim09_fixed[:,19],
      binned_sim10_fixed[:,19])
a20 = (binned_sim01_fixed[:,20],binned_sim02_fixed[:,20],binned_sim03_fixed[:,20],
      binned_sim04_fixed[:,20],binned_sim05_fixed[:,20],binned_sim06_fixed[:,20],
      binned_sim07_fixed[:,20],binned_sim08_fixed[:,20],binned_sim09_fixed[:,20],
      binned_sim10_fixed[:,20])
for i in range(10):
    plt.plot(a0[i], a19[i], 'b_', markersize = 3)
    plt.plot((a0[i],a0[i]), (a19[i]-a20[i],a19[i]+a20[i]), 'b', linewidth=.5, alpha = 0.2)
plt.plot(binned_e1dq1[:,0], binned_e1dq1[:,19], 'r_')
plt.savefig("./datasets/LC_3C454.3/dataplotfixed.png", dpi = 300, bbox_inches='tight')

# It is possible to do the simple mean of the variables because fluxes, times
# and backgruond are the same in all simulation and the only parameter
# interested in is Col.17 needed to find N_s
binned_mean = (binned_sim01_fixed+binned_sim02_fixed+binned_sim03_fixed+
 binned_sim04_fixed+binned_sim05_fixed+binned_sim06_fixed+
 binned_sim07_fixed+binned_sim08_fixed+binned_sim09_fixed+
 binned_sim10_fixed)/10
              
# Plot of binned_mean with its error and the real data
plt.plot(binned_e1dq1[:,0],binned_e1dq1[:,19], 'r_')
plt.plot(binned_mean[:,0],binned_mean[:,19], 'b_')
plt.plot((binned_mean[:,0],binned_mean[:,0]),
         (binned_mean[:,19]-binned_mean[:,20],binned_mean[:,19]+binned_mean[:,20]),
         'b', linewidth=.5, alpha = 0.2)
plt.savefig("./datasets/LC_3C454.3/RealDatavsMean.png", dpi = 300, bbox_inches='tight')


# Choose one of model
data = copy.copy(binned_mean)
N_b = data[:,26]
N_s = data[:,17] * data[:,2]
for i in range(len(N_s)):
    if N_s[i] < 0:
        N_s[i] = 0
np.mean(N_b),np.std(N_b)
np.mean(N_s),np.std(N_s)
print(stats.chisquare((data[:,19]-np.mean(data[:,19]))/np.std(data[:,19]),(binned_e1dq1[:,19]-np.mean(binned_e1dq1[:,19]))/np.std(binned_e1dq1[:,19])))
# SNR is a measure that compares the level of a desired signal to the level of background noise.
SNR = N_s/np.sqrt(N_s+2*N_b)
# S is for significance and it is an other way to evaluate relation between N_s and N:b
S = np.sqrt(2*((N_s+N_b)*np.log(2*(N_s+N_b)/(N_s+2*N_b))+N_b*np.log(2*(N_b)/(N_s+2*N_b))))
t = data[:,0] # tempo
x = np.round((N_s+N_b)/data[:,2]*1e+07) # frequenza eventi per ogni intervallo di tempo
N = len(x) # length of our dataset
p0 = 0.05 #false-positive rate
#Function of Bayesian block
edges = bayesian_blocks(t, x, fitness='events', gamma = np.exp(-(4-73.53*p0*N**(-0.478))))
#The output of the function is a vector of edges
edges
plt.hist(t, bins = edges)
plt.savefig("./datasets/LC_3C454.3/bbfixed.png", dpi = 300, bbox_inches='tight')
# It is important to work with the density inside the block
# a2[0] = Density in each block
# a2_cum is the cumulative density
a2 = np.histogram(t, bins = edges)
a2_cum = np.cumsum(a2[0])


#Define parameters within Bayesian block
c = 0 #counts
(b2,b3,b4,b5,b6) = (list(),list(),list(),list(),list()) 
# b2: vector of flux of energies grouped by blocks
# b3: N_s grouped by blocks
# b4: N_b grouped by blocks
# b5: exp grouped by blocks
# b6: error flux energy grouped by blocks
#Here we are going to build sets inside a vector related to each block
for j in range(len(a2_cum)):
    b2.append(data[c:a2_cum[j],19])
    b3.append(N_s[c:a2_cum[j]])
    b4.append(N_b[c:a2_cum[j]])
    b5.append(data[c:a2_cum[j],2])
    b6.append(data[c:a2_cum[j],20])
    c = a2_cum[j]
#Then we are create values needed to the plot creation
mean_block = list() # mean of flux within
prop_err_block = list() # the mean of the square root of the sum of col.20 squared
prop_upper_err_block = list() #mean+1sigma error
prop_lower_err_block = list() #mean-1sigma error
sum_exp_block = list() # sum of exposition inside the block
sum_N_b_block = list() # sum of N_b inside the block
sum_N_s_block = list() # sum of N_s inside the block
SNR_block = list() # SNR inside the block
rate_block = list() # rate inside the block
S_block = list() # Significance inside the block
cts_block = list()
for i in range(len(b2)):
    cts_block.append(np.mean(b3[i])+np.mean(b4[i]))
    S_block.append(np.sqrt(2)*np.sqrt((np.sum(b3[i])+np.sum(b4[i]))*np.log(2*(np.sum(b3[i])+np.sum(b4[i]))/(np.sum(b3[i])+2*np.sum(b4[i])))+np.sum(b4[i])*np.log(2*np.sum(b4[i])/(np.sum(b3[i])+2*np.sum(b4[i])))))
    rate_block.append((np.sum(b3[i])+np.sum(b4[i]))/np.sum(b5[i]))
    SNR_block.append(np.sum(b3[i])/np.sqrt(np.sum(b3[i])+np.sum(b4[i])))
    sum_N_s_block.append(sum(b3[i]))
    sum_N_b_block.append(sum(b4[i]))
    sum_exp_block.append(sum(b5[i]))
    mean_block.append(np.mean(b2[i]))
    prop_err_block.append(np.sqrt(np.sum(b6[i]**2))/a2[0][i])
    prop_upper_err_block.append(mean_block[i]+prop_err_block[i])
    prop_lower_err_block.append(mean_block[i]-prop_err_block[i])
    
# Here we are going to build a Moving average with lag equal to 8
sma = {'time' : data[:,0],'flux':data[:,19], 'var_flux':data[:,20]**2}
df = pd.DataFrame(sma)
df['SMA_8'] = df.iloc[:,1].rolling(window=8).mean()
df['SMA_8_err'] = np.sqrt(df.iloc[:,2].rolling(window=8).sum())/8


#### Plot Bayesian Block Above Data ####
# Bayesian blocks exept the edges (first and last)
plt.plot(((edges[1:]+edges[:-1])/2-(edges[1:]-edges[:-1])/2,(edges[1:]+edges[:-1])/2+(edges[1:]-edges[:-1])/2),(mean_block,mean_block), 'g', linewidth=1)
# First Bayesian block
plt.plot((edges[0],edges[0]-21600/2),(mean_block[0],mean_block[0]), 'g', linewidth=1)
# Last bayesian block
plt.plot((edges[len(edges)-1],edges[len(edges)-1]+21600/2),(mean_block[len(mean_block)-1],mean_block[len(mean_block)-1]), 'g', linewidth=1)
# Bayesian blocks' error
plt.plot(((edges[1:]+edges[:-1])/2,(edges[1:]+edges[:-1])/2),(prop_upper_err_block,prop_lower_err_block), 'g', linewidth=.5, alpha = 0.5)
# Data
plt.plot(data[:,0], data[:,19], 'b_', markersize = 3)
# Error of data (col.20)
plt.plot((data[:,0],data[:,0]), (data[:,19]-data[:,20],data[:,19]+data[:,20]), 'b', linewidth=.5, alpha = 0.5)
#Connection between the block in order to create a type of bars
plt.plot((edges[1:-1],edges[1:-1]),(mean_block[:-1],mean_block[1:]), 'g', linewidth=1)
plt.plot((edges[0]-21600/2,edges[0]-21600/2),(mean_block[0],0), 'g', linewidth=1)
plt.plot((edges[len(edges)-1]+21600/2,edges[len(edges)-1]+21600/2),(mean_block[-1],0), 'g', linewidth=1)
plt.xlabel('TT time *10$^8$')
plt.ylabel('Flux ph/cm2*s *10$^-5$')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#plt.savefig("./datasets/LC_3C454.3/plotfixed.png", dpi = 300, bbox_inches='tight')


#### Bayesian Block Above Data with Moving Average ####
# We are working with the moving average and the first 6 values of MA are nan, so
# in order to get a nice plot we remove the first 6 values of data and first 4 blocks

# Bayesian blocks exept the edges (first and last)
plt.plot(((edges[5:]+edges[4:-1])/2-(edges[5:]-edges[4:-1])/2,(edges[5:]+edges[4:-1])/2+(edges[5:]-edges[4:-1])/2),(mean_block[4:],mean_block[4:]), 'g', linewidth=1)
# First Bayesian block
plt.plot((edges[4],edges[4]-21600/2),(mean_block[4],mean_block[4]), 'g', linewidth=1)
# Last bayesian block
plt.plot((edges[len(edges)-1],edges[len(edges)-1]+21600/2),(mean_block[len(mean_block)-1],mean_block[len(mean_block)-1]), 'g', linewidth=1)
# Bayesian blocks' error
plt.plot(((edges[5:]+edges[4:-1])/2,(edges[5:]+edges[4:-1])/2),(prop_upper_err_block[4:],prop_lower_err_block[4:]), 'g', linewidth=.5, alpha = 0.5)
# Data
plt.plot(data[7:,0], data[7:,19], 'b_', markersize = 3)
# Error of data (col.20)
plt.plot((data[7:,0],data[7:,0]), (data[7:,19]-data[7:,20],data[7:,19]+data[7:,20]), 'b', linewidth=.5, alpha = 0.5)
#Connection between the block in order to create a type of bars
plt.plot((edges[5:-1],edges[5:-1]),(mean_block[4:-1],mean_block[5:]), 'g', linewidth=1)
plt.plot((edges[4]-21600/2,edges[4]-21600/2),(mean_block[4],min(data[:,19])), 'g', linewidth=1)
plt.plot((edges[len(edges)-1]+21600/2,edges[len(edges)-1]+21600/2),(mean_block[-1],min(data[:,19])), 'g', linewidth=1)
#MA(8) 2 days
plt.plot(df['time'][7:],df['SMA_8'][7:], 'r', linewidth = 1)
plt.plot(df['time'][7:],df['SMA_8'][7:]+3*df['SMA_8_err'][7:], 'r--', linewidth = 1)
# Background flux
plt.plot(df['time'][7:],data[7:,26]/data[7:,2])
# Condition that gives us which values are over the 3 sigma Confidence Interval
for i in range(len(data[:,0])):
    if data[i,19] > df['SMA_8'][i]+3*df['SMA_8_err'][i]:
        plt.plot(data[i,0], data[i,19], 'k_', markersize = 3)
        plt.plot((data[i,0],data[i,0]), (data[i,19]-data[i,20],data[i,19]+data[i,20]), 'k', linewidth=.5, alpha = 0.5)
plt.xlabel('TT time *10$^8$')
plt.ylabel('Flux ph/cm2*s *10$^-5$')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.savefig("./datasets/LC_3C454.3/smafixed.png", dpi = 300, bbox_inches='tight')


#### PLOT SNR ###
# Value of SNR
plt.plot(data[:,0], SNR, 'b_', markersize = 3)
# SNR error
plt.plot((data[:,0],data[:,0]), (SNR-np.sqrt(SNR),SNR+np.sqrt(SNR)), 'b', linewidth=.5, alpha = 0.5)
# Bayesian blocks of SNR (depends on the numbers of N_s and N_b inside the block)
plt.plot(((edges[1:]+edges[:-1])/2-(edges[1:]-edges[:-1])/2,(edges[1:]+edges[:-1])/2+(edges[1:]-edges[:-1])/2),(SNR_block,SNR_block), 'g', linewidth=1)
plt.plot((edges[0],edges[0]-21600/2),(SNR_block[0],SNR_block[0]), 'g', linewidth=1)
plt.plot((edges[len(edges)-1],edges[len(edges)-1]+21600/2),(SNR_block[len(SNR_block)-1],SNR_block[len(SNR_block)-1]), 'g', linewidth=1)
# Bayesian blocks' error
plt.plot(((edges[1:]+edges[:-1])/2,(edges[1:]+edges[:-1])/2),(SNR_block+np.sqrt(SNR_block),SNR_block-np.sqrt(SNR_block)), 'g', linewidth=.5, alpha = 0.5)
# Connection between the block in order to create a type of bars
plt.plot((edges[1:-1],edges[1:-1]),(SNR_block[:-1],SNR_block[1:]), 'g', linewidth=1)
plt.plot((edges[0]-21600/2,edges[0]-21600/2),(SNR_block[0],0), 'g', linewidth=1)
# Threshold
plt.plot(data[:,0],3*np.ones(len(data[:,0])), 'r:', linewidth = 1)
# Condition that gives us which values are over the  threshold = 3
for i in range(len(SNR_block)):
    if SNR_block[i] > 3:
        plt.plot(((edges[i+1]+edges[i])/2-(edges[i+1]-edges[i])/2,(edges[i+1]+edges[i])/2+(edges[i+1]-edges[i])/2),(SNR_block[i],SNR_block[i]), 'k', linewidth=1)
        plt.plot(((edges[i+1]+edges[i])/2,(edges[i+1]+edges[i])/2),(SNR_block[i]+np.sqrt(SNR_block[i]),SNR_block[i]-np.sqrt(SNR_block[i])), 'k', linewidth=.5, alpha = 0.5)
plt.title('Signal Noise Ratio')
plt.savefig("./datasets/LC_3C454.3/SNRblockfixed.png", dpi = 300, bbox_inches='tight')


#### PLOT Significance ####
# Value of S
plt.plot(data[:,0], S, 'b_', markersize = 3)
# S error
plt.plot((data[:,0],data[:,0]), (S-np.sqrt(S),S+np.sqrt(S)), 'b', linewidth=.5, alpha = 0.5)
# Bayesian blocks of S (depends on the numbers of N_s and N_b inside the block)
plt.plot(((edges[1:]+edges[:-1])/2-(edges[1:]-edges[:-1])/2,(edges[1:]+edges[:-1])/2+(edges[1:]-edges[:-1])/2),(S_block,S_block), 'g', linewidth=1)
plt.plot((edges[0],edges[0]-21600/2),(S_block[0],S_block[0]), 'g', linewidth=1)
plt.plot((edges[len(edges)-1],edges[len(edges)-1]+21600/2),(S_block[len(S_block)-1],S_block[len(S_block)-1]), 'g', linewidth=1)
# Bayesian blocks' error
plt.plot(((edges[1:]+edges[:-1])/2,(edges[1:]+edges[:-1])/2),(S_block+np.sqrt(S_block),S_block-np.sqrt(S_block)), 'g', linewidth=.5, alpha = 0.5)
# Connection between the block in order to create a type of bars
plt.plot((edges[1:-1],edges[1:-1]),(S_block[:-1],S_block[1:]), 'g', linewidth=1)
plt.plot((edges[0]-21600/2,edges[0]-21600/2),(S_block[0],0), 'g', linewidth=1)
# Threshold
plt.plot(data[:,0],3*np.ones(len(data[:,0])), 'r:', linewidth = 1)
# Condition that gives us which values are over the  threshold = 3
for i in range(len(S_block)):
    if S_block[i] > 3:
        plt.plot(((edges[i+1]+edges[i])/2-(edges[i+1]-edges[i])/2,(edges[i+1]+edges[i])/2+(edges[i+1]-edges[i])/2),(S_block[i],S_block[i]), 'k', linewidth=1)
        plt.plot(((edges[i+1]+edges[i])/2,(edges[i+1]+edges[i])/2),(S_block[i]+np.sqrt(S_block[i]),S_block[i]-np.sqrt(S_block[i])), 'k', linewidth=.5, alpha = 0.5)
plt.title('Significance')
plt.savefig("./datasets/LC_3C454.3/Sblockfixed.png", dpi = 300, bbox_inches='tight')

# The sum of exps in each block
plt.plot((edges[1:]+edges[:-1])/2,sum_exp_block) 
plt.savefig("./datasets/LC_3C454.3/expblockfixed.png", dpi = 300, bbox_inches='tight')

# Plot of signal flux and noise flux
plt.plot(data[:,0],(data[:,26]/data[:,2]))
plt.plot(data[:,0],(N_s/data[:,2]))
plt.legend(('NoiseFlux', "SignalFlux"))
plt.savefig("./datasets/LC_3C454.3/NSfluxfixed.png", dpi = 300, bbox_inches='tight')

# SNR without Bayesian blocks and three theshold levels
plt.plot(data[:,0],SNR, 'b.')
plt.plot(data[:,0],np.ones(len(data[:,0])), 'r')
plt.plot(data[:,0],2*np.ones(len(data[:,0])), 'r--')
plt.plot(data[:,0],3*np.ones(len(data[:,0])), 'r:')
plt.savefig("./datasets/LC_3C454.3/SNRfixed.png", dpi = 300, bbox_inches='tight')


#### Efficiency ####
print(stats.chisquare(data[:,3],binned_e1dq1[:,3]))
lista = list()
for i in range(len(mean_block)):
    lista += [cts_block[i]]*a2[0][i]
print(stats.chisquare(data[:,3],lista))

# >>> print(stats.chisquare(binned_mean[:,3],binned_e1dq1[:,3]))
# Power_divergenceResult(statistic=49.73015904855378, pvalue=0.9999999898064731)
# >>> lista = list()
# >>> for i in range(len(mean_block)):
# >>>     lista += [cts_block[i]]*a2[0][i]
# >>> print(stats.chisquare(data[:,3],lista))
# Power_divergenceResult(statistic=20.551938047737043, pvalue=1.0)

def RealTime(data):
    N_b = data[:,26]
    N_s = data[:,17] * data[:,2]
    #Needed to remove negative value of N_s
    for i in range(len(N_s)):
        if N_s[i] < 0:
            N_s[i] = 0
    # SNR is a measure that compares the level of a desired signal to the level of background noise.
    SNR = N_s/np.sqrt(N_s+2*N_b)
    # S is for significance and it is an other way to evaluate relation between N_s and N:b
    S = np.sqrt(2*((N_s+N_b)*np.log(2*(N_s+N_b)/(N_s+2*N_b))+N_b*np.log(2*(N_b)/(N_s+2*N_b))))
    # cumilative density within blocks
    a2_cum = list() 
    # Creation of moving average with 8 lags
    sma = {'time' : data[:,0],'flux':data[:,19], 'var_flux':data[:,20]**2}
    df = pd.DataFrame(sma)
    df['SMA_8'] = df.iloc[:,1].rolling(window=8).mean()
    df['SMA_8_err'] = np.sqrt(df.iloc[:,2].rolling(window=8).sum())/8
    # Begin the for cicle
    for i in range(len(data[:,0])):
        # Fix data in a scatter plot
        plt.plot(data[:i+1,0], data[:i+1,19], 'b_')
        y_03 = data[:i+1,19]-data[:i+1,20]
        y_04 = data[:i+1,19]+data[:i+1,20]
        plt.plot((data[:i+1,0][:i+1],data[:i+1,0][:i+1]),(y_03[:i+1],y_04[:i+1]), 'b', linewidth = 1, alpha = 0.3)
        plt.xlabel('TT time *10^8')
        plt.ylabel('Flux ph/cm2*s *10^-8') 
        if i > 1:
            # Bayesian block parameters
            t = data[:i+1,0] ## time
            x = np.round((N_s[:i+1]+N_b[:i+1])/data[:i+1,2]*1e+07) ## event
            # Bayesian block function that gives as output edges of blocks
            edges = bayesian_blocks(t, x, fitness='events', gamma = np.exp(-(4-73.53*0.05*(i+1)**(-0.478)))) #bayesian blocks function
            # a2[0] gives us the density in each block
            # a2[1] is also the edges of blocks
            a2 = np.histogram(t, bins = edges)
            a2_cum = np.cumsum(a2[0])
            # some parameter needed to plot our data
            # c = count
            # b2: vector of flux of energies grouped by blocks
            # b3: error flux energy grouped by blocks
            (c,b2,b3) = (0,list(),list())
            # mean_block: mean of flux within
            # prop_err_block: the mean of the square root of the sum of col.20 squared
            # prop_upper_err_block: mean+1sigma error
            # prop_lower_err_block: mean-1sigma error
            (mean_block,prop_err_block,prop_upper_err_block,prop_lower_err_block) = (list(),list(),list(),list())
            for j in range(len(a2_cum)):
                b2.append(data[c:a2_cum[j],19]) # data[left_edge_block:right_edge_block] 
                b3.append(data[c:a2_cum[j],20])
                mean_block.append(np.mean(b2[j]))
                prop_err_block.append(np.sqrt(np.sum(b3[j]**2))/a2[0][j])
                prop_upper_err_block.append(mean_block[j]+prop_err_block[j])
                prop_lower_err_block.append(mean_block[j]-prop_err_block[j])  
                # Bayesian block using as height the mean_block
                plt.plot((t[c]-21600/2,t[a2_cum[j]-1]+21600/2),(mean_block[-1],mean_block[-1]), 'g', linewidth=1) ## value
                # Bayesian blocks' error
                plt.plot(((t[c]+t[a2_cum[j]-1])/2,(t[c]+t[a2_cum[j]-1])/2),(prop_upper_err_block[j],prop_lower_err_block[j]), 'g', linewidth=1, alpha = 0.3)
                c = a2_cum[j]
            if len(mean_block)>1:
                # Connection between the block in order to create a type of bars
                plt.plot((edges[1:-1],edges[1:-1]),(mean_block[:-1],mean_block[1:]), 'g', linewidth=1)
                plt.plot((t[0]-21600/2,t[0]-21600/2),(0,mean_block[0]), 'g', linewidth=1)
                plt.plot((t[a2_cum[-1]-1]+21600/2,t[a2_cum[-1]-1]+21600/2),(0,mean_block[-1]), 'g', linewidth=1)
        # Red lines are the MA(8) and the +3sigma CI
        plt.plot(df['time'][7:i+1],df['SMA_8'][7:i+1], 'r', linewidth = 1)
        plt.plot(df['time'][7:i+1],df['SMA_8'][7:i+1]+3*df['SMA_8_err'][7:i+1], 'r--', linewidth = 1)
        plt.xlabel('TT time *10$^8$')
        plt.ylabel('Flux ph/cm2*s *10$^-5$')
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.show()
        if i > 1:
            (c,b4,b5,SNR_block) = (0,list(),list(),list())
            for j in range(len(a2_cum)):
                b4.append(N_b[c:a2_cum[j]])  
                b5.append(N_s[c:a2_cum[j]])
                SNR_block.append(np.sum(b5[j])/np.sqrt(np.sum(b5[j])+np.sum(b4[j])))
                c = a2_cum[j]
            if len(SNR_block)>1:
                # Value of SNR
                plt.plot(data[:i+1,0], SNR[:i+1], 'b_', markersize = 3)
                # SNR error
                plt.plot((data[:i+1,0],data[:i+1,0]), (SNR[:i+1]-np.sqrt(SNR[:i+1]),SNR[:i+1]+np.sqrt(SNR[:i+1])), 'b', linewidth=.5, alpha = 0.5)
                # Bayesian blocks of SNR (depends on the numbers of N_s and N_b inside the block)
                plt.plot(((edges[1:]+edges[:-1])/2-(edges[1:]-edges[:-1])/2,(edges[1:]+edges[:-1])/2+(edges[1:]-edges[:-1])/2),(SNR_block,SNR_block), 'g', linewidth=1)
                plt.plot((edges[0],edges[0]-21600/2),(SNR_block[0],SNR_block[0]), 'g', linewidth=1)
                plt.plot((edges[len(edges)-1],edges[len(edges)-1]+21600/2),(SNR_block[len(SNR_block)-1],SNR_block[len(SNR_block)-1]), 'g', linewidth=1)
                # Bayesian blocks' error
                plt.plot(((edges[1:]+edges[:-1])/2,(edges[1:]+edges[:-1])/2),(SNR_block+np.sqrt(SNR_block),SNR_block-np.sqrt(SNR_block)), 'g', linewidth=.5, alpha = 0.5)
                # Connection between the block in order to create a type of bars
                plt.plot((edges[1:-1],edges[1:-1]),(SNR_block[:-1],SNR_block[1:]), 'g', linewidth=1)
                plt.plot((edges[0]-21600/2,edges[0]-21600/2),(SNR_block[0],0), 'g', linewidth=1)
                plt.plot((edges[-1]+21600/2,edges[-1]+21600/2),(0,SNR_block[-1]), 'g', linewidth=1)
                # Threshold
                plt.plot(data[:i+1,0],3*np.ones(len(data[:i+1,0])), 'r:', linewidth = 1)
                # Condition that gives us which values are over the  threshold = 3
                for ij in range(len(SNR_block)):
                    if SNR_block[ij] > 3:
                        plt.plot(((edges[ij+1]+edges[ij])/2-(edges[ij+1]-edges[ij])/2,(edges[ij+1]+edges[ij])/2+(edges[ij+1]-edges[ij])/2),(SNR_block[ij],SNR_block[ij]), 'k', linewidth=1)
                        plt.plot(((edges[ij+1]+edges[ij])/2,(edges[ij+1]+edges[ij])/2),(SNR_block[ij]+np.sqrt(SNR_block[ij]),SNR_block[ij]-np.sqrt(SNR_block[ij])), 'k', linewidth=.5, alpha = 0.5)
                plt.title('Signal Noise Ratio')
                plt.xlabel('TT time *10$^8$')
                plt.ylabel('SNR')
                plt.show()
            (c,b4,b5,S_block) = (0,list(),list(),list())
            for j in range(len(a2_cum)):
                b4.append(N_b[c:a2_cum[j]])  
                b5.append(N_s[c:a2_cum[j]])
                S_block.append(np.sqrt(2)*np.sqrt((np.sum(b5[j])+np.sum(b4[j]))*np.log(2*(np.sum(b5[j])+np.sum(b4[j]))/(np.sum(b5[j])+2*np.sum(b4[j])))+np.sum(b4[j])*np.log(2*np.sum(b4[j])/(np.sum(b5[j])+2*np.sum(b4[j])))))
                c = a2_cum[j]
            if len(S_block)>1:
                # Value of S
                plt.plot(data[:i+1,0], S[:i+1], 'b_', markersize = 3)
                # S error
                plt.plot((data[:i+1,0],data[:i+1,0]), (S[:i+1]-np.sqrt(S[:i+1]),S[:i+1]+np.sqrt(S[:i+1])), 'b', linewidth=.5, alpha = 0.5)
                # Bayesian blocks of S (depends on the numbers of N_s and N_b inside the block)
                plt.plot(((edges[1:]+edges[:-1])/2-(edges[1:]-edges[:-1])/2,(edges[1:]+edges[:-1])/2+(edges[1:]-edges[:-1])/2),(S_block,S_block), 'g', linewidth=1)
                plt.plot((edges[0],edges[0]-21600/2),(S_block[0],S_block[0]), 'g', linewidth=1)
                plt.plot((edges[len(edges)-1],edges[len(edges)-1]+21600/2),(S_block[len(S_block)-1],S_block[len(S_block)-1]), 'g', linewidth=1)
                # Bayesian blocks' error
                plt.plot(((edges[1:]+edges[:-1])/2,(edges[1:]+edges[:-1])/2),(S_block+np.sqrt(S_block),S_block-np.sqrt(S_block)), 'g', linewidth=.5, alpha = 0.5)
                # Connection between the block in order to create a type of bars
                plt.plot((edges[1:-1],edges[1:-1]),(S_block[:-1],S_block[1:]), 'g', linewidth=1)
                plt.plot((edges[0]-21600/2,edges[0]-21600/2),(S_block[0],0), 'g', linewidth=1)
                plt.plot((edges[-1]+21600/2,edges[-1]+21600/2),(0,S_block[-1]), 'g', linewidth=1)
                # Threshold
                plt.plot(data[:i+1,0],3*np.ones(len(data[:i+1,0])), 'r:', linewidth = 1)
                # Condition that gives us which values are over the  threshold = 3
                for ij in range(len(S_block)):
                    if S_block[ij] > 3:
                        plt.plot(((edges[ij+1]+edges[ij])/2-(edges[ij+1]-edges[ij])/2,(edges[ij+1]+edges[ij])/2+(edges[ij+1]-edges[ij])/2),(S_block[ij],S_block[ij]), 'k', linewidth=1)
                        plt.plot(((edges[ij+1]+edges[ij])/2,(edges[ij+1]+edges[ij])/2),(S_block[ij]+np.sqrt(S_block[ij]),S_block[ij]-np.sqrt(S_block[ij])), 'k', linewidth=.5, alpha = 0.5)
                plt.title('S')
                plt.xlabel('TT time *10$^8$')
                plt.ylabel('S')
                plt.show()
        if i < 116:
            print('Do you want to stop the program? Press y to quit')
            stop = input()
            if stop == 'y' or stop == 'Y':
                break