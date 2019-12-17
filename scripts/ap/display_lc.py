import os,sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from matplotlib.gridspec import GridSpec


def find_min_max_axes(array_one,array_two):

        min_one = np.amin(array_one)
        min_two = np.amin(array_two)
        min = np.amin([min_two,min_one])

        max_one = np.amax(array_one)
        max_two = np.amax(array_two)
        max = np.amax([max_one,max_two])

        return min,max


def extract_data(file_name):

        detection_array = []

        with open(file_name) as fp:
           for cnt, line in enumerate(fp):
               sqrtts = 0
               if(file_name.endswith(".lc")):

                   tstart = float(line.split()[13])
                   tstop =  float(line.split()[14])
                   flux =  float(line.split()[0])*flux_notation
                   exp = float(line.split()[8])/exp_notation
                   exp_norm = -1
                   flux_err =  float(line.split()[1])*flux_notation
                   if(float(line.split()[7])>0):
                       sqrtts =  float(line.split()[7])
                   else:
                       sqrtts =  -1
                   count = -1
                   count_bkg = -1
                   count_err = -1
                   count_bkg_err = -1
                   rate = -1
                   rate_err = -1

               if(file_name.endswith(".ap3")):

                   tstart = float(line.split()[0])
                   tstop =  float(line.split()[1])
                   rate = float(line.split()[19])*flux_notation
                   if(rate < 0.0):
                       print(line)
                       #continue
                   rate_err = float(line.split()[20])*flux_notation
                   if(rate_err > 2000.0):
                       print(line)
                       #continue
                       
                   flux =  float(line.split()[21])*flux_notation
                   if(flux < 0.0):
                       print(line)
                       #continue
                   exp = float(line.split()[2])/exp_notation
                   exp_norm = float(line.split()[7])/exp_notation
                   flux_err =  float(line.split()[22])*flux_notation
                   if(float(line.split()[23])>0):
                       sqrtts =  np.sqrt(float(line.split()[23]))
                   else:
                       sqrtts =  -1
                   count = float(line.split()[3])
                   count_err = float(np.sqrt(count))/2
                   count_bkg = float(line.split()[26])
                   count_bkg_err = float(np.sqrt(count_bkg))/2

               #FILTRO
#               if(sqrtts < 3.0):
#                   continue
#               else:
#                   print(file_name)
#                   print(line)

               detection_array.append({"rate":rate,"rate_err":rate_err,"count":count,"count_err":count_err,"count_bkg":count_bkg,"count_bkg_err":count_bkg_err,"tstart":tstart,"tstop":tstop,"flux":flux,"flux_err":flux_err,"sqrtts":sqrtts,"exp":exp,"exp_norm":exp_norm})


        return detection_array


# add filters in this funcion
def get_value_from_array(data_array):

    tstart_array = []
    tstop_array = []
    flux_array = []
    flux_err_array = []
    sqrtts_array = []
    x_array = []
    xerr_array = []
    exp_array = []
    exp_norm_array = []
    count_array = []
    count_bkg_array=[]
    count_err_array = []
    count_bkg_err_array=[]
    rate_array = []
    rate_err_array=[]

    for detection in data_array:

        tstart = detection['tstart']
        tstop = detection['tstop']
        flux = detection['flux']
        flux_err = detection['flux_err']
        sqrtts = detection['sqrtts']
        exp = detection['exp']
        exp_norm = detection['exp_norm']
        count = detection['count']
        count_bkg = detection['count_bkg']
        count_err = detection['count_err']
        count_bkg_err = detection['count_bkg_err']
        rate = detection['rate']
        rate_err = detection['rate_err']


        tstart_array.append(tstart)
        tstop_array.append(tstop)
        flux_array.append(flux)
        flux_err_array.append(flux_err)
        sqrtts_array.append(sqrtts)
        count_array.append(count)
        count_bkg_array.append(count_bkg)
        x_array.append(tstart+(tstop-tstart)/2)
        xerr_array.append((tstop-tstart)/2)
        exp_array.append(exp)
        exp_norm_array.append(exp_norm)
        count_err_array.append(count_err)
        count_bkg_err_array.append(count_bkg_err)
        rate_array.append(rate)
        rate_err_array.append(rate_err)




    return {"x":x_array,"xerr":xerr_array,"tstart":tstart_array,"tstop":tstop_array,"flux":flux_array,"flux_err":flux_err_array,"sqrtts":sqrtts_array,'exp':exp_array,'exp_norm':exp_norm_array,"count":count_array,'count_err':count_err_array,"count_bkg":count_bkg_array,'count_bkg_err':count_bkg_err_array,'rate':rate_array,'rate_err':rate_err_array}



###### MAIN #########

mode = sys.argv[1]
scatter = sys.argv[2]
file_one = sys.argv[3]
if(mode=="2"):
    file_two = sys.argv[4]

flux_notation = 100000000
exp_notation = 1000000


dict_one = extract_data(file_one)

binsize = 15

if(mode=="1"):


    data_array_one = get_value_from_array(dict_one)

    plt.figure(1)

    plt.errorbar(data_array_one['x'], data_array_one['flux'],xerr=data_array_one['xerr'], yerr=data_array_one['flux_err'], label='flux 1', fmt='r.')

    plt.xlabel("TT time *10^8")
    plt.ylabel("Flux ph/cm2 s *10^-8")
    plt.grid(True)
    plt.legend()

    plt.figure(1)
    fig, ax = plt.subplots(2,5)
    fig.suptitle('File 1 (red): '+os.path.basename(file_one))
    n, bins, patches = ax[0,0].hist(data_array_one['flux'], binsize, density=True, facecolor='r', alpha=0.6,label='flux 1')
    ax[0,0].plot(bins, norm.pdf(bins, np.mean(data_array_one['flux']), np.std(data_array_one['flux'])), color="black", linestyle="--", alpha=0.9)
    n, bins, patches = ax[0,1].hist(data_array_one['flux_err'], binsize, density=True, facecolor='r', alpha=0.6,label='flux err 1')
    ax[0,1].plot(bins, norm.pdf(bins, np.mean(data_array_one['flux_err']), np.std(data_array_one['flux_err'])), color="black", linestyle="--", alpha=0.9)
    n, bins, patches = ax[0,2].hist(data_array_one['sqrtts'], binsize, density=True, facecolor='r', alpha=0.6,label='sqrtts 1')
    ax[0,2].plot(bins, norm.pdf(bins, np.mean(data_array_one['sqrtts']), np.std(data_array_one['sqrtts'])), color="black", linestyle="--", alpha=0.9)
    ax[0,0].set_title("Flux")
    ax[0,0].set(xlabel='Flux ph/cm2 s *10^-8', ylabel='Det. Num.')
    ax[0,1].set_title("Flux err")
    ax[0,1].set(xlabel='Flux ph/cm2 s *10^-8')
    ax[0,2].set_title("sqrt(Ts)")
    ax[0,2].set(xlabel='sqrt(TS)')
    #plt.legend([l1, l2, l3],["HHZ 1", "HHN", "HHE"])

    ax[0,3].scatter(data_array_one['sqrtts'],data_array_one['flux'],color='red',s=10)
    ax[0,3].set(xlabel='sqrts(TS)', ylabel='Flux ph/cm2 s *10^-8')
    ax[0,3].set_title("sqrt(TS) vs Flux")

    n, bins, patches = ax[1,3].hist(data_array_one['exp'], binsize, density=True, facecolor='r', alpha=0.6,label='exp 1')
    ax[1,3].plot(bins, norm.pdf(bins, np.mean(data_array_one['exp']), np.std(data_array_one['exp'])), color="black", linestyle="--", alpha=0.9)
    ax[1,3].set_title("Exp")
    ax[1,3].set(xlabel='Exp 10^6')

    if(file_one.endswith(".ap3")):
        ax[0,4].scatter(data_array_one['exp_norm'],data_array_one['flux'],color='red',s=10)
        ax[0,4].set(xlabel='Exp Norm 10^6', ylabel='Flux ph/cm2 s *10^-8')
        ax[0,4].set_title("Exp norm vs Flux")

        n, bins, patches = ax[1,4].hist(data_array_one['exp_norm'], binsize, density=True, facecolor='r', alpha=0.6,label='exp 1')
        ax[1,4].plot(bins, norm.pdf(bins, np.mean(data_array_one['exp_norm']), np.std(data_array_one['exp_norm'])), color="black", linestyle="--", alpha=0.9)

        ax[1,4].set_title("Exp Norm")
        ax[1,4].set(xlabel='Exp Norm 10^6')

    plt.figure(1)
    fig, ax = plt.subplots(3,3)
    fig.suptitle('File 1 (red): '+os.path.basename(file_one))

    bins = 40

    #LC COUNT
    ax[0,0].errorbar(data_array_one['x'], data_array_one['count'],xerr=data_array_one['xerr'], yerr=data_array_one['count_err'], label=os.path.basename(file_one), fmt='r.')
    ax[0,0].set(xlabel="TT time *10^8", ylabel="Counts")
    ax[0,0].grid(True)

    #LC COUNT BKG
    ax[0,1].errorbar(data_array_one['x'], data_array_one['count_bkg'],xerr=data_array_one['xerr'], yerr=data_array_one['count_bkg_err'], label=os.path.basename(file_one), fmt='r.')
    ax[0,1].set(xlabel="TT time *10^8", ylabel="Counts Bkg")
    ax[0,1].grid(True)

    #LC EXP
    ax[1,0].errorbar(data_array_one['x'], data_array_one['exp'],xerr=data_array_one['xerr'],  label=os.path.basename(file_one), fmt='r.')
    ax[1,0].set(xlabel="TT time *10^8", ylabel="Exp")
    ax[1,0].grid(True)

    #LC EXP NORM
    ax[1,1].errorbar(data_array_one['x'], data_array_one['exp_norm'],xerr=data_array_one['xerr'] , label=os.path.basename(file_one), fmt='r.')
    ax[1,1].set(xlabel="TT time *10^8", ylabel="Exp Norm")
    ax[1,1].grid(True)


    #LC RATE
    ax[0,2].errorbar(data_array_one['x'], data_array_one['rate'],xerr=data_array_one['xerr'],yerr=data_array_one['rate_err'],  label=os.path.basename(file_one), fmt='r.')
    ax[0,2].set(xlabel="TT time *10^8", ylabel="Rate")
    ax[0,2].grid(True)

    #ISTOGRAMMA COUNT

    mean = np.mean(data_array_one['count'])
    std = np.std(data_array_one['count'])
    n, bins, patches = ax[2,2].hist(data_array_one['count'], binsize, density=True, facecolor='r', alpha=0.6,label="mean: "+str(round(mean,2))+",std= "+str(round(std,2)))
    ax[2,2].plot(bins, norm.pdf(bins, mean, std), color="black", linestyle="--", alpha=0.9)
    ax[2,2].set_title("Count")
    ax[2,2].set(xlabel='Count')
    ax[2,2].legend()



    #ISTOGRAMMA COUNT BKG
    mean = np.mean(data_array_one['count_bkg'])
    std = np.std(data_array_one['count_bkg'])
    n, bins, patches = ax[1,2].hist(data_array_one['count_bkg'], binsize, density=True, facecolor='r', alpha=0.6,label="mean: "+str(round(mean,2))+",std= "+str(round(std,2)))
    ax[1,2].plot(bins, norm.pdf(bins, mean, std), color="black", linestyle="--", alpha=0.9)
    ax[1,2].set_title("Count")
    ax[1,2].set(xlabel='Count')
    ax[1,2].legend()


    plt.show()

if(mode=="2"):

    dict_two = extract_data(file_two)

    print(len(dict_one))
    print(len(dict_two))

    if(scatter=="1"):
        plot_scatter = True
    else:
        plot_scatter = False


    if(plot_scatter is True):
        dict_two_clean = []
        dict_one_clean = []


        for detection_two in dict_two:
            tstart_two = detection_two['tstart']

            for detection_one in dict_one:
                tstart_one = detection_one['tstart']

                if(tstart_one == tstart_two):
                    dict_two_clean.append(detection_two)
                    dict_one_clean.append(detection_one)


        print(len(dict_one_clean))
        print(len(dict_two_clean))


        data_array_one = get_value_from_array(dict_one_clean)
        data_array_two = get_value_from_array(dict_two_clean)
    else:
        data_array_one = get_value_from_array(dict_one)
        data_array_two = get_value_from_array(dict_two)


    plt.figure(1)

    plt.errorbar(data_array_one['x'], data_array_one['flux'],xerr=data_array_one['xerr'], yerr=data_array_one['flux_err'], label=os.path.basename(file_one), fmt='r.')
    plt.errorbar(data_array_two['x'], data_array_two['flux'],xerr=data_array_two['xerr'], yerr=data_array_two['flux_err'], label=os.path.basename(file_two), fmt='b.')

    plt.xlabel("TT time *10^8")
    plt.ylabel("Flux ph/cm2 s *10^-8")
    plt.grid(True)
    plt.legend()

    plt.figure(1)
    fig, ax = plt.subplots(2,5)
    fig.suptitle('File 1 (red): '+os.path.basename(file_one)+', File 2 (blue): '+os.path.basename(file_two))
    n, bins, patches = ax[0,0].hist(data_array_one['flux'], binsize, density=True, facecolor='r', alpha=0.6,label='flux 1')
    ax[0,0].plot(bins, norm.pdf(bins, np.mean(data_array_one['flux']), np.std(data_array_one['flux'])), color="black", linestyle="--", alpha=0.9)
    n, bins, patches = ax[0,0].hist(data_array_two['flux'], binsize, density=True, facecolor='b', alpha=0.6,label='flux 2')
    ax[0,0].plot(bins, norm.pdf(bins, np.mean(data_array_two['flux']), np.std(data_array_two['flux'])), color="black", linestyle="--", alpha=0.9)

    n, bins, patches = ax[0,1].hist(data_array_one['flux_err'], binsize, density=True, facecolor='r', alpha=0.6,label='flux err 1')
    ax[0,1].plot(bins, norm.pdf(bins, np.mean(data_array_one['flux_err']), np.std(data_array_one['flux_err'])), color="black", linestyle="--", alpha=0.9)
    n, bins, patches = ax[0,1].hist(data_array_two['flux_err'], binsize, density=True, facecolor='b', alpha=0.6,label='flux err 2')
    ax[0,1].plot(bins, norm.pdf(bins, np.mean(data_array_two['flux_err']), np.std(data_array_two['flux_err'])), color="black", linestyle="--", alpha=0.9)

    n, bins, patches = ax[0,2].hist(data_array_one['sqrtts'], binsize, density=True, facecolor='r', alpha=0.6,label='sqrt(TS) 1')
    ax[0,2].plot(bins, norm.pdf(bins, np.mean(data_array_one['sqrtts']), np.std(data_array_one['sqrtts'])), color="black", linestyle="--", alpha=0.9)
    n, bins, patches = ax[0,2].hist(data_array_two['sqrtts'], binsize, density=True, facecolor='b', alpha=0.6,label='sqrt(TS) 2')
    ax[0,2].plot(bins, norm.pdf(bins, np.mean(data_array_two['sqrtts']), np.std(data_array_two['sqrtts'])), color="black", linestyle="--", alpha=0.9)

    ax[0,0].set_title("Flux")
    ax[0,0].set(xlabel='Flux ph/cm2 s *10^-8', ylabel='Det. Num.')
    ax[0,1].set_title("Flux err")
    ax[0,1].set(xlabel='Flux ph/cm2 s *10^-8')
    ax[0,2].set_title("sqrt(TS)")
    ax[0,2].set(xlabel='sqrt(TS)')


    ##### CHANGE SCATTER PLOT ##############################

    #change axes and label

    #possible axes: sqrtts,flux,flux_err,exp,exp_norm

    axes_x = "exp"
    x_label = "exp"

    axes_y = "exp_norm"
    y_label = "exp_norm"

    ax[0,3].scatter(data_array_one[axes_x],data_array_one[axes_y],color='red',s=10)
    ax[0,3].scatter(data_array_two[axes_x],data_array_two[axes_y],color='blue',s=10)
    ax[0,3].set(xlabel=x_label, ylabel=y_label)
    ax[0,3].set_title(axes_x+" vs "+axes_y)

    ###########################################################

    if(file_one.endswith(".ap3")):

        ax[0,4].scatter(data_array_one['exp_norm'],data_array_one['flux'],color='red',s=10)
        ax[0,4].set(xlabel='Exp Norm 10^6', ylabel='Flux ph/cm2 s *10^-8')
        ax[0,4].set_title("Exp norm vs Flux")

        n, bins, patches = ax[1,4].hist(data_array_one['exp_norm'], binsize, density=True, facecolor='r', alpha=0.6,label='exp 1')
        ax[1,4].plot(bins, norm.pdf(bins, np.mean(data_array_one['exp_norm']), np.std(data_array_one['exp_norm'])), color="black", linestyle="--", alpha=0.9)
        ax[1,4].set_title("Exp Norm")
        ax[1,4].set(xlabel='Exp Norm 10^6')

    if(file_two.endswith(".ap3")):

        ax[0,4].scatter(data_array_two['exp_norm'],data_array_two['flux'],color='blue',s=10)
        ax[0,4].set(xlabel='Exp Norm 10^6', ylabel='Flux ph/cm2 s *10^-8')
        ax[0,4].set_title("Exp Norm vs Flux")

        n, bins, patches = ax[1,4].hist(data_array_two['exp_norm'], binsize, density=True, facecolor='b', alpha=0.6,label='exp 1')
        ax[1,4].plot(bins, norm.pdf(bins, np.mean(data_array_two['exp_norm']), np.std(data_array_two['exp_norm'])), color="black", linestyle="--", alpha=0.9)
        ax[1,4].set_title("Exp Norm")
        ax[1,4].set(xlabel='Exp Norm 10^6')

    mean = np.mean(data_array_one['exp'])
    std = np.std(data_array_one['exp'])
    n, bins, patches = ax[1,3].hist(data_array_one['exp'], binsize, density=True, facecolor='r', alpha=0.6,label="mean: "+str(round(mean,2))+",std= "+str(round(std,2)))
    ax[1,3].plot(bins, norm.pdf(bins, mean, std), color="black", linestyle="--", alpha=0.9)

    mean = np.mean(data_array_two['exp'])
    std = np.std(data_array_two['exp'])
    n, bins, patches = ax[1,3].hist(data_array_two['exp'], binsize, density=True, facecolor='b', alpha=0.6,label="mean: "+str(round(mean,2))+",std= "+str(round(std,2)))
    ax[1,3].plot(bins, norm.pdf(bins, mean, std), color="black", linestyle="--", alpha=0.9)
    ax[1,3].set_title("Exp")
    ax[1,3].set(xlabel='Exp 10^6')
    ax[1,3].legend()



    if(plot_scatter):
        ax[1,0].scatter(data_array_one['flux'], data_array_two['flux'],s=10)
        min,max = find_min_max_axes(data_array_one['flux'],data_array_two['flux'])
        ax[1,0].set_xlim([min,max])
        ax[1,0].set_ylim([min,max])
        ax[1,0].plot([min,max],[min,max],color='red')

    if(plot_scatter):
        ax[1,1].scatter(data_array_one['flux_err'],data_array_two['flux_err'],s=10)
        min,max = find_min_max_axes(data_array_one['flux_err'],data_array_two['flux_err'])
        ax[1,1].set_xlim([min,max])
        ax[1,1].set_ylim([min,max])
        ax[1,1].plot([min,max],[min,max],color='red')

    if(plot_scatter):
        ax[1,2].scatter(data_array_one['sqrtts'],data_array_two['sqrtts'],s=10)
        min,max = find_min_max_axes(data_array_one['sqrtts'],data_array_two['sqrtts'])
        ax[1,2].set_xlim([min,max])
        ax[1,2].set_ylim([min,max])
        ax[1,2].plot([min,max],[min,max],color='red')

    #ax[1,0].set_title("Flux")
    ax[1,0].set(xlabel='Flux 1', ylabel='Flux 2')
    #ax[1,1].set_title("Flux err")
    ax[1,1].set(xlabel='Flux err 1',ylabel="Flux err 2")
    #ax[1,2].set_title("Ts")
    ax[1,2].set(xlabel='sqrt(TS) 1',ylabel='sqrt(TS) 2')


    plt.figure(1)
    fig, ax = plt.subplots(3,3)
    fig.suptitle('File 1 (red): '+os.path.basename(file_one)+', File 2 (blue): '+os.path.basename(file_two))

    bins = 40

    #LC COUNT
    ax[0,0].errorbar(data_array_one['x'], data_array_one['count'],xerr=data_array_one['xerr'], yerr=data_array_one['count_err'], label=os.path.basename(file_one), fmt='r.')
    ax[0,0].errorbar(data_array_two['x'], data_array_two['count'],xerr=data_array_two['xerr'], yerr=data_array_two['count_err'], label=os.path.basename(file_two), fmt='b.')
    ax[0,0].set(xlabel="TT time *10^8", ylabel="Counts")
    ax[0,0].grid(True)

    #LC COUNT BKG
    ax[0,1].errorbar(data_array_one['x'], data_array_one['count_bkg'],xerr=data_array_one['xerr'], yerr=data_array_one['count_bkg_err'], label=os.path.basename(file_one), fmt='r.')
    ax[0,1].errorbar(data_array_two['x'], data_array_two['count_bkg'],xerr=data_array_two['xerr'], yerr=data_array_two['count_bkg_err'], label=os.path.basename(file_two), fmt='b.')
    ax[0,1].set(xlabel="TT time *10^8", ylabel="Counts Bkg")
    ax[0,1].grid(True)

    #LC EXP
    ax[1,0].errorbar(data_array_one['x'], data_array_one['exp'],xerr=data_array_one['xerr'],  label=os.path.basename(file_one), fmt='r.')
    ax[1,0].errorbar(data_array_two['x'], data_array_two['exp'],xerr=data_array_two['xerr'],  label=os.path.basename(file_two), fmt='b.')
    ax[1,0].set(xlabel="TT time *10^8", ylabel="Exp")
    ax[1,0].grid(True)

    #LC EXP NORM
    ax[1,1].errorbar(data_array_one['x'], data_array_one['exp_norm'],xerr=data_array_one['xerr'] , label=os.path.basename(file_one), fmt='r.')
    ax[1,1].errorbar(data_array_two['x'], data_array_two['exp_norm'],xerr=data_array_two['xerr'], label=os.path.basename(file_two), fmt='b.')
    ax[1,1].set(xlabel="TT time *10^8", ylabel="Exp Norm")
    ax[1,1].grid(True)

    #### SCATER PLOT COUNT
    if(plot_scatter):
        ax[2,0].scatter(data_array_one['count'],data_array_two['count'],s=10)
        min,max = find_min_max_axes(data_array_one['count'],data_array_two['count'])
        ax[2,0].set_xlim([min,max])
        ax[2,0].set_ylim([min,max])
        ax[2,0].plot([min,max],[min,max],color='red')
        ax[2,0].set_title("Count vs Count")

    #### SCATER PLOT COUNT BKG
    if(plot_scatter):
        ax[2,1].scatter(data_array_one['count_bkg'],data_array_two['count_bkg'],s=10)
        min,max = find_min_max_axes(data_array_one['count_bkg'],data_array_two['count_bkg'])
        ax[2,1].set_xlim([min,max])
        ax[2,1].set_ylim([min,max])
        ax[2,1].plot([min,max],[min,max],color='red')
        ax[2,1].set_title("Count bkg vs Count bkg")



    #LC RATE
    ax[0,2].errorbar(data_array_one['x'], data_array_one['rate'],xerr=data_array_one['xerr'],yerr=data_array_one['rate_err'],  label=os.path.basename(file_one), fmt='r.')
    ax[0,2].errorbar(data_array_two['x'], data_array_two['rate'],xerr=data_array_two['xerr'],yerr=data_array_two['rate_err'] ,  label=os.path.basename(file_two), fmt='b.')
    ax[0,2].set(xlabel="TT time *10^8", ylabel="Rate")
    ax[0,2].grid(True)

    #ISTOGRAMMA COUNT

    mean = np.mean(data_array_one['count'])
    std = np.std(data_array_one['count'])
    n, bins, patches = ax[2,2].hist(data_array_one['count'], binsize, density=True, facecolor='r', alpha=0.6,label="mean: "+str(round(mean,2))+",std= "+str(round(std,2)))
    ax[2,2].plot(bins, norm.pdf(bins, mean, std), color="black", linestyle="--", alpha=0.9)
    ax[2,2].set_title("Count")
    ax[2,2].set(xlabel='Count')
    ax[2,2].legend()


    mean = np.mean(data_array_two['count'])
    std = np.std(data_array_two['count'])
    n, bins, patches = ax[2,2].hist(data_array_two['count'], binsize, density=True, facecolor='b', alpha=0.6,label="mean: "+str(round(mean,2))+",std= "+str(round(std,2)))
    ax[2,2].plot(bins, norm.pdf(bins, mean, std), color="black", linestyle="--", alpha=0.9)
    ax[2,2].set_title("Count")
    ax[2,2].set(xlabel='Count')
    ax[2,2].legend()


    #ISTOGRAMMA COUNT BKG
    mean = np.mean(data_array_one['count_bkg'])
    std = np.std(data_array_one['count_bkg'])
    n, bins, patches = ax[1,2].hist(data_array_one['count_bkg'], binsize, density=True, facecolor='r', alpha=0.6,label="mean: "+str(round(mean,2))+",std= "+str(round(std,2)))
    ax[1,2].plot(bins, norm.pdf(bins, mean, std), color="black", linestyle="--", alpha=0.9)
    ax[1,2].set_title("Count")
    ax[1,2].set(xlabel='Count')
    ax[1,2].legend()


    mean = np.mean(data_array_two['count_bkg'])
    std = np.std(data_array_two['count_bkg'])
    n, bins, patches = ax[1,2].hist(data_array_two['count_bkg'], binsize, density=True, facecolor='b', alpha=0.6,label="mean: "+str(round(mean,2))+",std= "+str(round(std,2)))
    ax[1,2].plot(bins, norm.pdf(bins, mean, std), color="black", linestyle="--", alpha=0.9)
    ax[1,2].set_title("Count")
    ax[1,2].set(xlabel='Count')
    ax[1,2].legend()

    plt.show()
