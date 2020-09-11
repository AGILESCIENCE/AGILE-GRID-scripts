# DESCRIPTION
#       Agileap: AGILE Observatory Aperture Photometry Analysis
# NOTICE
#      Any information contained in this software
#      is property of the AGILE TEAM and is strictly
#      private and confidential.
#      Copyright (C) 2005-2020 AGILE Team.
#          Bulgarelli Andrea <andrea.bulgarelli@inaf.it>
#          Guerrini Marco <marcoguerrini1993@gmail.com>
#          Antonio Addis <antonio.addis@inaf.it>
#          Valentina Fioretti <valentina.fioretti@inaf.it>
#          Parmiggiani Nicolò <nicolo.parmiggiani@inaf.it>
#      All rights reserved.
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
from astropy.stats import bayesian_blocks
import matplotlib.pyplot as plt
import pandas as pd
from tkinter import Tk
from tkinter.filedialog import askopenfilename
import copy
from astropy.io import fits 
from decimal import Decimal 

class BayesianBlocks:

    def BBdataset(self, data):
        t = data[:, 0]  # tempo
        esponente = '%.2E' % Decimal(min(data[:, 1]))
        if esponente[-3:] == '-01':
            # frequenza eventi per ogni intervallo di tempo
            x = np.round(data[:, 1]*1e+02)
        else:
            # frequenza eventi per ogni intervallo di tempo
            x = np.round(data[:, 1]*1e+07)
        #x = np.round(data[:,1]/np.abs(min(data[:,1])))
        N = len(x)  # length of our dataset
        p0 = 0.05  # false-positive rate
        #Function of Bayesian block
        edges = bayesian_blocks(t, x, fitness='events',
                                gamma=np.exp(-(4-73.53*p0*N**(-0.478))))
        a2 = np.histogram(t, bins=edges)
        a2_cum = np.cumsum(a2[0])
        # Define and create parameters within Bayesian block
        c = 0  # counts
        (b2, b3, mean_block, err_block, lista1, lista2) = (
            list(), list(), list(), list(), list(), list())
        # b2: vector of flux of energies grouped by blocks
        # b3: error flux energy grouped by blocks
        # mean of flux within
        # the mean of the square root of the sum of col.20 squared
        # lista1 = vector of height (len=117)
        # lista2 = vector of error (len=117)
        for j in range(len(a2_cum)):
            b2.append(data[c:a2_cum[j], 1])
            b3.append(data[c:a2_cum[j], 2])
            mean_block.append(np.mean(b2[j]))
            err_block.append(np.sqrt(np.sum(b3[j]**2))/a2[0][j])
            lista1 += [mean_block[j]]*a2[0][j]
            lista2 += [err_block[j]]*a2[0][j]
            c = a2_cum[j]
        m = edges[:-1]  # tstart
        m = np.append(m, edges[1:])  # tstop
        m = np.append(m, mean_block)  # height
        m = np.append(m, err_block)  # error
        m = np.transpose(m.reshape(4, len(mean_block)))
        m = pd.DataFrame(m)
        m = m.rename(columns={0: "tstart", 1: "tstop",
                            2: "FluxMeanBB", 3: "FluxErrorBB"})
        m.to_csv('BBdataset.csv', index=True, sep=';')
        d = data[:, 0]  # tstart
        d = np.append(d, data[:, 1])  # tstop
        d = np.append(d, lista1)  # height
        d = np.append(d, lista2)  # error
        d = np.transpose(d.reshape(4, len(lista1)))
        d = pd.DataFrame(d)
        d = d.rename(columns={0: "tstart", 1: "tstop",
                            2: "FluxMeanBB", 3: "FluxErrorBB"})
        d.to_csv('BBdataset2.csv', index=True, sep=';')
        return(m)

    def AllData(self, data, N_b, N_s):
        # SNR is a measure that compares the level of a desired signal to the level of background noise.
        # Formula of SNR = N_s/np.sqrt(N_s+2*N_b)
        # S is for significance and it is an other way to evaluate relation between N_s and N:b
        # Formula of S = np.sqrt(2*((N_s+N_b)*np.log(2*(N_s+N_b)/(N_s+2*N_b))+N_b*np.log(2*(N_b)/(N_s+2*N_b))))
        # cumilative density within blocks
        t = data[:, 0]  # tempo
        #x = np.round(data[:,1]/np.abs(min(data[:,1])))
        # frequenza eventi per ogni intervallo di tempo
        x = np.round(data[:, 1]*1e+07)
    #Function of Bayesian block
        edges = bayesian_blocks(t, x, fitness='events', gamma=3*1e-07)
        # It is important to work with the density inside the block
        # a2[0] = Density in each block
        # a2_cum is the cumulative density
        a2 = np.histogram(t, bins=edges)
        a2_cum = np.cumsum(a2[0])
        c = 0  # counts
        (b2, b3, b4, b5, b6) = (list(), list(), list(), list(), list())
        # b2: vector of flux of energies grouped by blocks
        # b3: N_s grouped by blocks
        # b4: N_b grouped by blocks
        # b5: exp grouped by blocks
        # b6: error flux energy grouped by blocks
        #Here we are going to build sets inside a vector related to each block
        #plt.ion()
        for j in range(len(a2_cum)):
            b2.append(data[c:a2_cum[j], 1])
            b3.append(N_s[c:a2_cum[j]])
            b4.append(N_b[c:a2_cum[j]])
            b5.append(data[c:a2_cum[j], 3])
            b6.append(data[c:a2_cum[j], 2])
            c = a2_cum[j]
        #Then we are create values needed to the plot creation
        mean_block = list()  # mean of flux within
        prop_err_block = list()  # the mean of the square root of the sum of col.20 squared
        prop_upper_err_block = list()  # mean+1sigma error
        prop_lower_err_block = list()  # mean-1sigma error
        sum_N_b_block = list()  # sum of N_b inside the block
        sum_N_s_block = list()  # sum of N_s inside the block
        SNR_block = list()  # SNR inside the block
        S_block = list()  # Significance inside the block
        fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(12, 16))
        for i in range(len(b2)):
            S_block.append((np.sqrt(2)*np.sqrt((np.sum(b3[i])+np.sum(b4[i]))*np.log(2*(np.sum(b3[i])+np.sum(b4[i]))/(np.sum(
                b3[i])+2*np.sum(b4[i])))+np.sum(b4[i])*np.log(2*np.sum(b4[i])/(np.sum(b3[i])+2*np.sum(b4[i])))))/sum(b5[i]))
            SNR_block.append(
                (np.sum(b3[i])/np.sqrt(np.sum(b3[i])+2*np.sum(b4[i])))/sum(b5[i]))
            sum_N_s_block.append(sum(b3[i]))
            sum_N_b_block.append(sum(b4[i]))
            mean_block.append(np.mean(b2[i]))
            prop_err_block.append(np.sqrt(np.sum(b6[i]**2))/a2[0][i])
            prop_upper_err_block.append(mean_block[i]+prop_err_block[i])
            prop_lower_err_block.append(mean_block[i]-prop_err_block[i])
        # Bayesian blocks exept the edges (first and last)
        ax0.plot(((edges[1:]+edges[:-1])/2-(edges[1:]-edges[:-1])/2, (edges[1:]+edges[:-1]
                                                                    )/2+(edges[1:]-edges[:-1])/2), (mean_block, mean_block), 'g', linewidth=1)
        # First Bayesian block
        ax0.plot((edges[0], edges[0]-21600/2),
                (mean_block[0], mean_block[0]), 'g', linewidth=1)
        # Last bayesian block
        ax0.plot((edges[len(edges)-1], edges[len(edges)-1]+21600/2),
                (mean_block[len(mean_block)-1], mean_block[len(mean_block)-1]), 'g', linewidth=1)
        # Bayesian blocks' error
        ax0.plot(((edges[1:]+edges[:-1])/2, (edges[1:]+edges[:-1])/2),
                (prop_upper_err_block, prop_lower_err_block), 'g', linewidth=.5, alpha=0.5)
        # Data
        ax0.plot(data[:, 0], data[:, 1], 'b_', markersize=3)
        # Error of data (col.20)
        ax0.plot((data[:, 0], data[:, 0]), (data[:, 1]-data[:, 2],
                                            data[:, 1]+data[:, 2]), 'b', linewidth=.5, alpha=0.5)
        #Connection between the block in order to create a type of bars
        ax0.plot((edges[1:-1], edges[1:-1]),
                (mean_block[:-1], mean_block[1:]), 'g', linewidth=1)
        ax0.plot((edges[0]-21600/2, edges[0]-21600/2),
                (mean_block[0], 0), 'g', linewidth=1)
        ax0.plot((edges[len(edges)-1]+21600/2, edges[len(edges)-1] +
                21600/2), (mean_block[-1], 0), 'g', linewidth=1)
        ax0.set_xlabel('TT time *10$^8$')
        ax0.set_ylabel('Flux ph/cm2*s *10$^-5$')
        ax0.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        #plt.show()
        #SNR
        # Bayesian blocks of SNR (depends on the numbers of N_s and N_b inside the block)
        ax1.plot(((edges[1:]+edges[:-1])/2-(edges[1:]-edges[:-1])/2, (edges[1:]+edges[:-1]
                                                                    )/2+(edges[1:]-edges[:-1])/2), (SNR_block, SNR_block), 'g', linewidth=1)
        ax1.plot((edges[0], edges[0]-21600/2),
                (SNR_block[0], SNR_block[0]), 'g', linewidth=1)
        ax1.plot((edges[len(edges)-1], edges[len(edges)-1]+21600/2),
                (SNR_block[len(SNR_block)-1], SNR_block[len(SNR_block)-1]), 'g', linewidth=1)
        # Connection between the block in order to create a type of bars
        ax1.plot((edges[1:-1], edges[1:-1]),
                (SNR_block[:-1], SNR_block[1:]), 'g', linewidth=1)
        ax1.plot((edges[0]-21600/2, edges[0]-21600/2),
                (SNR_block[0], 0), 'g', linewidth=1)
        ax1.plot((edges[len(edges)-1]+21600/2, edges[len(edges)-1] +
                21600/2), (SNR_block[len(SNR_block)-1], 0), 'g', linewidth=1)
        # Threshold
        ax1.plot(data[:, 0], (np.mean(N_s)/np.mean(N_s+2*N_b)) /
                np.mean(data[:, 3])*5*np.ones(len(data)), 'r:', linewidth=1)
        # Condition that gives us which values are over the  threshold = 3
        count = 0
        for i in range(len(SNR_block)):
            if SNR_block[i] > (np.mean(N_s)/np.mean(N_s+2*N_b))/np.mean(data[:, 3])*5:
                ax1.plot(((edges[i+1]+edges[i])/2-(edges[i+1]-edges[i])/2, (edges[i+1]+edges[i])/2+(
                    edges[i+1]-edges[i])/2), (SNR_block[i], SNR_block[i]), 'k', linewidth=1)
                count += 1
        ax1.set_title('Signal Noise Ratio (SNR)')
        ax1.set_xlabel('TT time *10$^8$')
        ax1.set_ylabel('SNR')
        print(count)
        #plt.show()
        #### PLOT Significance S ####
        # Bayesian blocks of S (depends on the numbers of N_s and N_b inside the block)
        ax2.plot(((edges[1:]+edges[:-1])/2-(edges[1:]-edges[:-1])/2, (edges[1:] +
                                                                    edges[:-1])/2+(edges[1:]-edges[:-1])/2), (S_block, S_block), 'g', linewidth=1)
        ax2.plot((edges[0], edges[0]-21600/2),
                (S_block[0], S_block[0]), 'g', linewidth=1)
        ax2.plot((edges[len(edges)-1], edges[len(edges)-1]+21600/2),
                (S_block[len(S_block)-1], S_block[len(S_block)-1]), 'g', linewidth=1)
        # Bayesian blocks' error
        # Connection between the block in order to create a type of bars
        ax2.plot((edges[1:-1], edges[1:-1]),
                (S_block[:-1], S_block[1:]), 'g', linewidth=1)
        ax2.plot((edges[0]-21600/2, edges[0]-21600/2),
                (S_block[0], 0), 'g', linewidth=1)
        ax2.plot((edges[len(edges)-1]+21600/2, edges[len(edges)-1] +
                21600/2), (S_block[len(S_block)-1], 0), 'g', linewidth=1)
        # Threshold
        ax2.plot(data[:, 0], (np.mean(N_s)/np.mean(N_s+2*N_b)) /
                np.mean(data[:, 3])*5*np.ones(len(data)), 'r:', linewidth=1)
        # Condition that gives us which values are over the  threshold = 3
        count = 0
        for i in range(len(S_block)):
            if S_block[i] > (np.mean(N_s)/np.mean(N_s+2*N_b))/np.mean(data[:, 3])*5:
                count += 1
                ax2.plot(((edges[i+1]+edges[i])/2-(edges[i+1]-edges[i])/2, (edges[i+1]+edges[i])/2+(
                    edges[i+1]-edges[i])/2), (S_block[i], S_block[i]), 'k', linewidth=1)
        ax2.set_title('Significance (S)')
        ax2.set_xlabel('TT time *10$^8$')
        ax2.set_ylabel('S')
        print(count)
        plt.show()

    def RealTime(self, data, N_b, N_s):
        # SNR is a measure that compares the level of a desired signal to the level of background noise.
        # SNR = N_s/np.sqrt(N_s+2*N_b)
        # S is for significance and it is an other way to evaluate relation between N_s and N:b
        # S = np.sqrt(2*((N_s+N_b)*np.log(2*(N_s+N_b)/(N_s+2*N_b))+N_b*np.log(2*(N_b)/(N_s+2*N_b))))
        # cumilative density within blocks
        for i in range(len(data[:, 0])):
            # Fix data in a scatter plot
            plt.plot(data[:i+1, 0], data[:i+1, 1], 'b_')
            y_03 = data[:i+1, 1]-data[:i+1, 2]
            y_04 = data[:i+1, 1]+data[:i+1, 2]
            plt.plot((data[:i+1, 0][:i+1], data[:i+1, 0][:i+1]),
                    (y_03[:i+1], y_04[:i+1]), 'b', linewidth=1, alpha=0.3)
            plt.xlabel('TT time *10^8')
            plt.ylabel('Flux ph/cm2*s *10^-8')
            if i > 1:
                # Bayesian block parameters
                t = data[:i+1, 0]  # time
                #x = np.round(data[:i+1,1]/np.abs(min(data[:i+1,1])))
                x = np.round(data[:i+1, 1]*1e+07)  # event
                # Bayesian block function that gives as output edges of blocks
                # bayesian blocks function
                edges = bayesian_blocks(t, x, fitness='events', gamma=3*1e-07)
                # a2[0] gives us the density in each block
                # a2[1] is also the edges of blocks
                a2 = np.histogram(t, bins=edges)
                a2_cum = np.cumsum(a2[0])
                # some parameter needed to plot our data
                # c = count
                # b2: vector of flux of energies grouped by blocks
                # b6: error flux energy grouped by blocks
                (c, b2, b6) = (0, list(), list())
                # mean_block: mean of flux within
                # prop_err_block: the mean of the square root of the sum of col.20 squared
                # prop_upper_err_block: mean+1sigma error
                # prop_lower_err_block: mean-1sigma error
                (mean_block, prop_err_block, prop_upper_err_block,
                prop_lower_err_block) = (list(), list(), list(), list())
                for j in range(len(a2_cum)):
                    # data[left_edge_block:right_edge_block]
                    b2.append(data[c:a2_cum[j], 1])
                    b6.append(data[c:a2_cum[j], 2])
                    mean_block.append(np.mean(b2[j]))
                    prop_err_block.append(np.sqrt(np.sum(b6[j]**2))/a2[0][j])
                    prop_upper_err_block.append(mean_block[j]+prop_err_block[j])
                    prop_lower_err_block.append(mean_block[j]-prop_err_block[j])
                    # Bayesian block using as height the mean_block
                    plt.plot((t[c]-21600/2, t[a2_cum[j]-1]+21600/2),
                            (mean_block[-1], mean_block[-1]), 'g', linewidth=1)  # value
                    # Bayesian blocks' error
                    plt.plot(((t[c]+t[a2_cum[j]-1])/2, (t[c]+t[a2_cum[j]-1])/2),
                            (prop_upper_err_block[j], prop_lower_err_block[j]), 'g', linewidth=1, alpha=0.3)
                    c = a2_cum[j]
                if len(mean_block) > 1:
                    # Connection between the block in order to create a type of bars
                    plt.plot((edges[1:-1], edges[1:-1]),
                            (mean_block[:-1], mean_block[1:]), 'g', linewidth=1)
                    plt.plot((t[0]-21600/2, t[0]-21600/2),
                            (0, mean_block[0]), 'g', linewidth=1)
                    plt.plot((t[a2_cum[-1]-1]+21600/2, t[a2_cum[-1]-1] +
                            21600/2), (0, mean_block[-1]), 'g', linewidth=1)
            plt.xlabel('TT time *10$^8$')
            plt.ylabel('Flux ph/cm2*s *10$^-5$')
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            plt.show()
            if i > 1:
                (c, b3, b4, b5, SNR_block) = (0, list(), list(), list(), list())
                # b3: N_s grouped by blocks
                # b4: N_b grouped by blocks
                # b5: exp grouped by blocks
                for j in range(len(a2_cum)):
                    b4.append(N_b[c:a2_cum[j]])
                    b3.append(N_s[c:a2_cum[j]])
                    b5.append(data[c:a2_cum[j], 3])
                    SNR_block.append(
                        (np.sum(b3[j])/np.sqrt(np.sum(b3[j])+2*np.sum(b4[j])))/sum(b5[j]))
                    c = a2_cum[j]
                if len(SNR_block) > 1:
                    # Plot of SNR
                    plt.plot(((edges[1:]+edges[:-1])/2-(edges[1:]-edges[:-1])/2, (edges[1:]+edges[:-1])/2+(
                        edges[1:]-edges[:-1])/2), (SNR_block, SNR_block), 'g', linewidth=1)
                    plt.plot((edges[0], edges[0]-21600/2),
                            (SNR_block[0], SNR_block[0]), 'g', linewidth=1)
                    plt.plot((edges[len(edges)-1], edges[len(edges)-1]+21600/2),
                            (SNR_block[len(SNR_block)-1], SNR_block[len(SNR_block)-1]), 'g', linewidth=1)
                    # Connection between the block in order to create a type of bars
                    plt.plot((edges[1:-1], edges[1:-1]),
                            (SNR_block[:-1], SNR_block[1:]), 'g', linewidth=1)
                    plt.plot((edges[0]-21600/2, edges[0]-21600/2),
                            (SNR_block[0], 0), 'g', linewidth=1)
                    plt.plot((edges[-1]+21600/2, edges[-1]+21600/2),
                            (0, SNR_block[-1]), 'g', linewidth=1)
                    # Threshold
                    plt.plot(data[:i+1, 0], (np.mean(N_s[:i+1])/np.mean(N_s[:i+1]+2*N_b[:i+1])) /
                            np.mean(data[:i+1, 3])*5*np.ones(len(data[:i+1, 0])), 'r:', linewidth=1)
                    # Condition that gives us which values are over the  threshold = 3
                    for ij in range(len(SNR_block)):
                        if SNR_block[ij] > (np.mean(N_s[:i+1])/np.mean(N_s[:i+1]+2*N_b[:i+1]))/np.mean(data[:i+1, 3])*5:
                            plt.plot(((edges[ij+1]+edges[ij])/2-(edges[ij+1]-edges[ij])/2, (edges[ij+1]+edges[ij])/2+(
                                edges[ij+1]-edges[ij])/2), (SNR_block[ij], SNR_block[ij]), 'k', linewidth=1)
                    plt.title('Signal Noise Ratio')
                    plt.xlabel('TT time *10$^8$')
                    plt.ylabel('SNR')
                    plt.show()
                (c, b3, b4, b5, S_block) = (0, list(), list(), list(), list())
                for j in range(len(a2_cum)):
                    b4.append(N_b[c:a2_cum[j]])
                    b3.append(N_s[c:a2_cum[j]])
                    b5.append(data[c:a2_cum[j], 3])
                    S_block.append(np.sqrt(2)*np.sqrt((np.sum(b3[j])+np.sum(b4[j]))*np.log(2*(np.sum(b3[j])+np.sum(b4[j]))/(np.sum(
                        b3[j])+2*np.sum(b4[j])))+np.sum(b4[j])*np.log(2*np.sum(b4[j])/(np.sum(b3[j])+2*np.sum(b4[j]))))/sum(b5[j]))
                    c = a2_cum[j]
                if len(S_block) > 1:
                    # Value of S
                    plt.plot(((edges[1:]+edges[:-1])/2-(edges[1:]-edges[:-1])/2, (edges[1:]+edges[:-1])/2+(
                        edges[1:]-edges[:-1])/2), (S_block, S_block), 'g', linewidth=1)
                    plt.plot((edges[0], edges[0]-21600/2),
                            (S_block[0], S_block[0]), 'g', linewidth=1)
                    plt.plot((edges[len(edges)-1], edges[len(edges)-1]+21600/2),
                            (S_block[len(S_block)-1], S_block[len(S_block)-1]), 'g', linewidth=1)
                    # Connection between the block in order to create a type of bars
                    plt.plot((edges[1:-1], edges[1:-1]),
                            (S_block[:-1], S_block[1:]), 'g', linewidth=1)
                    plt.plot((edges[0]-21600/2, edges[0]-21600/2),
                            (S_block[0], 0), 'g', linewidth=1)
                    plt.plot((edges[-1]+21600/2, edges[-1]+21600/2),
                            (0, S_block[-1]), 'g', linewidth=1)
                    # Threshold
                    plt.plot(data[:i+1, 0], (np.mean(N_s[:i+1])/np.mean(N_s[:i+1]+2*N_b[:i+1])) /
                            np.mean(data[:i+1, 3])*5*np.ones(len(data[:i+1, 0])), 'r:', linewidth=1)
                    # Condition that gives us which values are over the  threshold = 3
                    for ij in range(len(S_block)):
                        if S_block[ij] > (np.mean(N_s[:i+1])/np.mean(N_s[:i+1]+2*N_b[:i+1]))/np.mean(data[:i+1, 3])*5:
                            plt.plot(((edges[ij+1]+edges[ij])/2-(edges[ij+1]-edges[ij])/2, (edges[ij+1]+edges[ij])/2+(
                                edges[ij+1]-edges[ij])/2), (S_block[ij], S_block[ij]), 'k', linewidth=1)
                    plt.title('S')
                    plt.xlabel('TT time *10$^8$')
                    plt.ylabel('S')
                    plt.show()
            if i < 116:
                print('Do you want to stop the program? Press y to quit')
                stop = input()
                if stop == 'y' or stop == 'Y':
                    break
    
    
    def main(self, filepath, sel=2, values=100):
        a = filepath

        ###---AP3
        if a[len(a)-4:len(a)] == '.ap3':
            dataset = np.fromfile(a, sep=' ')
            dataset = dataset.reshape(len(dataset)//27, 27)
            dataset = dataset[2:, ]
            data = copy.copy(dataset)
            N_b = data[:, 26]
            N_s = data[:, 17] * data[:, 2]
            for i in range(len(N_s)):
                if N_s[i] < 0:
                    N_s[i] = 0
            m = data[:, 0]
            m = np.append(m, (N_s+N_b)/data[:, 2])
            m = np.append(m, data[:, 20])
            m = np.append(m, data[:, 2])
            m = np.transpose(m.reshape(4, 117))
        
        ###----CSV
        elif a[len(a)-4:len(a)] == '.csv':  # da correggere
            data_input = pd.read_csv(a, sep=';')
            dataset = data_input.values
            dataset = dataset[:, 1:]
            if len(dataset[1]) == 27:
                data = copy.copy(dataset)
                N_b = data[:, 26]
                N_s = data[:, 17] * data[:, 2]
                m = data[:, 0]
                m = np.append(m, (N_s+N_b)/data[:, 2])
                m = np.append(m, data[:, 20])
                m = np.append(m, data[:, 2])
                m = np.transpose(m.reshape(4, 117))
        
        ###----FITS
        elif a[len(a)-4:len(a)] == 'fits':  # da correggere
            hdulist = fits.open(a)
            hdu_data = hdulist[1].data
            data_tsorted = np.sort(hdu_data, order='TIME')
            t_tsorted = data_tsorted.field(1)
            work1 = t_tsorted
            h = np.histogram(work1, range=[0, np.int(
                max(work1))+1], bins=np.int(max(work1))+1)
            N_tot = h[0]
            rate = (N_tot)/(5*5*np.pi)  # exp fissa
            rate_err = np.sqrt(rate)
            m = range(np.int(max(work1))+1)
            m = np.append(m, rate)
            m = np.append(m, rate_err)
            m = np.transpose(m.reshape(3, len(range(np.int(max(work1))+1))))
            self.BBdataset(m)
        
        ###---AP4
        elif a[len(a)-4:len(a)] == '.ap4':
            dataset = pd.read_csv(a, sep=" ", header=0)
            N_b = dataset["cts_rateWeightedMeanR4"]
            N_s = dataset["ratediffR4"] * dataset["exp"]
            N_s.loc[N_s < 0] = 0
            col1 = (N_s + N_b)/dataset["exp"]
            m = dataset[["tstart", "rateError", "exp"]]
            m.insert(1, "col1", col1, True)
            m = m.to_numpy()
            m = m[:values, :]

        if sel == 1:
            self.BBdataset(m)
        elif sel == 2:
            self.AllData(m, N_b, N_s)
        elif sel == 3:
            self.RealTime(m, N_b, N_s)
