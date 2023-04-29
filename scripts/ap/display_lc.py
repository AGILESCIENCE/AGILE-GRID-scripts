import os,sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator

#### USAGE #####

# python display_lc.py [mode 1-2] [scatter plot 0-1] [fixed flux or -1] [tstart_window_mjd] [tstart_window_mjd] filepath1 filepath2 [optional, only if mode=2]
# mode 1 -> only one file, mode = 2 -> two file in input
# fixed flux showed inside plot


def time_tt_to_mjd(tt_time):
     return (float(tt_time) / 86400.0)+53005.0

def find_min_max_axes(array_one,array_two):

        min_one = np.amin(array_one)
        min_two = np.amin(array_two)
        min = np.amin([min_two,min_one])

        max_one = np.amax(array_one)
        max_two = np.amax(array_two)
        max = np.amax([max_one,max_two])

        return min,max


def extract_data(file_name,tstart_window_mjd,tstop_window_mjd):

        detection_array = []

        with open(file_name) as fp:
            for cnt, line in enumerate(fp):
                sqrtts = 0

                if(file_name.endswith(".fermi.lc")):
                    #folder sources ts npred flux flux_err flux_ul95 flux100_ul95 eflux eflux_err eflux100 eflux100_err eflux_ul95 sens3 sens35 sens4 sens5 tmin tmax tmin(MJD) tmax(MJD) phase
                    tstart = float(line.split()[19])
                    tstop =  float(line.split()[20])

                    #check time window
                    if(tstart_window_mjd!=-1 and tstop_window_mjd!=-1):
                       if(tstart<tstart_window_mjd or tstart>tstop_window_mjd ):
                           continue


                    x = tstart+(tstop-tstart)/2
                    x_err =  (tstop-tstart)/2
                    flux =  float(line.split()[4])*flux_notation
                    if(flux == -1):
                       flux=0
                    exp = 0
                    exp_norm = -1
                    flux_err =  float(line.split()[5])*flux_notation
                    if(float(line.split()[2])>0):
                       sqrtts =  np.sqrt(float(line.split()[2]))
                    else:
                       sqrtts =  0
                    count = -1
                    count_bkg = -1
                    count_err = -1
                    count_bkg_err = -1
                    rate = -1
                    rate_err = -1
                    flux_ul = float(line.split()[6])*flux_notation
                    sensitivity = -1

                    if(sqrtts<3):
                       flux_ul = flux

                if(file_name.endswith(".mle.lc")):

                    tstart = float(line.split()[11])
                    tstop =  float(line.split()[12])

                    #check time window
                    if(tstart_window_mjd!=-1 and tstop_window_mjd!=-1):
                       if(tstart<tstart_window_mjd or tstart>tstop_window_mjd ):
                           continue


                    x = tstart+(tstop-tstart)/2
                    x_err =  (tstop-tstart)/2
                    flux =  float(line.split()[0])*flux_notation
                    if(flux == -1):
                       flux=0
                    exp = float(line.split()[8])/exp_notation
                    exp_norm = -1
                    flux_err =  float(line.split()[1])*flux_notation
                    if(float(line.split()[7])>0):
                       sqrtts =  float(line.split()[7])
                    else:
                       sqrtts =  0
                    count = -1
                    count_bkg = -1
                    count_err = -1
                    count_bkg_err = -1
                    rate = -1
                    rate_err = -1
                    flux_ul = -1
                    sensitivity = -1

                    if(sqrtts<3):
                       flux_ul = flux

                if(file_name.endswith(".ap3") or file_name.endswith(".ap4")):

                    if(line.startswith("tstart")):
                       continue

                    components = line.split()

                    tstart = time_tt_to_mjd(float(components[0]))
                    tstop =  time_tt_to_mjd(float(components[1]))

                    print(tstart)

                    #check time window
                    if(tstart_window_mjd!=-1 and tstop_window_mjd!=-1):
                        if(tstart<tstart_window_mjd or tstart>tstop_window_mjd ):
                            print("continue")
                            continue


                    x = tstart+(tstop-tstart)/2
                    x_err =  (tstop-tstart)/2
                    rate = float(components[19])*flux_notation
                    #if(rate < 0.0):
                    #    print(line)
                       #continue
                    rate_err = float(components[20])*flux_notation
                    #if(rate_err >   2000.0):
                    #    print(line)
                       #continue


                    flux =  float(components[21])*flux_notation
                    #if(flux < 0.0):
                    #    print(line)
                       #continue
                    exp = float(components[2])/exp_notation
                    exp_norm = float(components[7])/exp_notation
                    flux_err =  float(components[22])*flux_notation
                    #23 Sa
                    #27 Slima
                    sindex=27
                    if(float(components[sindex])>0):
                       sqrtts =  float(components[sindex]) #S
                    else:
                       sqrtts =  0
                    count = float(components[3])
                    count_err = float(np.sqrt(count))
                    count_bkg = float(components[26])
                    count_bkg_err = float(np.sqrt(count_bkg))

                    flux_ul = float(components[31])*flux_notation
                    #31 -> 2
                    #35 -> 3
                    #39 -> 4
                    sensitivity = float(components[39])*flux_notation


                    #FILTRO
                    #if(sqrtts < 3):
                    #  flux_err = 0
                     # flux=flux_ul
                    #  continue
                    #else:
                    #  print(file_name)
                      #print(line)

                detection_array.append({"x":x,"x_err":x_err,"rate":rate,"rate_err":rate_err,"count":count,"count_err":count_err,"count_bkg":count_bkg,"count_bkg_err":count_bkg_err,"tstart":tstart,"tstop":tstop,"flux":flux,"flux_err":flux_err,"sqrtts":sqrtts,"exp":exp,"exp_norm":exp_norm,"flux_ul":flux_ul,"sensitivity":sensitivity})


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
    flux_ul_array = []
    fixed_array = []

    for detection in data_array:

        tstart = detection['tstart']
        tstop = detection['tstop']
        x = detection['x']
        x_err = detection['x_err']
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
        flux_ul = detection["flux_ul"]
        sensitivity = detection["sensitivity"]


        tstart_array.append(tstart)
        tstop_array.append(tstop)
        flux_array.append(flux)
        flux_err_array.append(flux_err)
        sqrtts_array.append(sqrtts)
        count_array.append(count)
        count_bkg_array.append(count_bkg)
        x_array.append(x)
        xerr_array.append(x_err)
        exp_array.append(exp)
        exp_norm_array.append(exp_norm)
        count_err_array.append(count_err)
        count_bkg_err_array.append(count_bkg_err)
        rate_array.append(rate)
        rate_err_array.append(rate_err)
        flux_ul_array.append(flux_ul)
        fixed_array.append(sensitivity)


    return {"x":x_array,"xerr":xerr_array,"tstart":tstart_array,"tstop":tstop_array,"flux":flux_array,"flux_err":flux_err_array,"sqrtts":sqrtts_array,'exp':exp_array,'exp_norm':exp_norm_array,"count":count_array,'count_err':count_err_array,"count_bkg":count_bkg_array,'count_bkg_err':count_bkg_err_array,'rate':rate_array,'rate_err':rate_err_array,'sensitivity':fixed_array,'flux_ul':flux_ul_array}


def main():
    ###### MAIN #########

    mode = sys.argv[1]
    scatter = sys.argv[2]
    fixed_flux = float(sys.argv[3])
    tstart_window_mjd = float(sys.argv[4])
    tstop_window_mjd = float(sys.argv[5])

    file_one = sys.argv[6]
    if(mode=="2"):
        file_two = sys.argv[7]

    flux_notation = 100000000
    exp_notation = 1000000


    dict_one = extract_data(file_one,tstart_window_mjd,tstop_window_mjd)

    binsize = 15

    if(mode=="1"):


        data_array_one = get_value_from_array(dict_one)


        fig0 = plt.figure(constrained_layout=False)

        count = 0
        for detection in dict_one:
            if(count==0):
                count=1
                continue
            
            if(detection['sqrtts']<3):
                #add arrow
                #plt.arrow(detection['x'],detection['flux'],0,-70,width=1,head_width=10,head_starts_at_zero=True)
                plt.errorbar(detection['x'], detection['flux_ul'],xerr=detection['x_err'], yerr=0, fmt='rv')
            else:
                plt.errorbar(detection['x'], detection['flux'],xerr=detection['x_err'], yerr=detection['flux_err'], fmt='r.')

        if(dict_one[0]['sqrtts']<3):
            plt.errorbar(dict_one[0]['x'], dict_one[0]['flux_ul'],xerr=dict_one[0]['x_err'], yerr=0, label=os.path.basename(file_one), fmt='rv')
        else:
            plt.errorbar(dict_one[0]['x'], dict_one[0]['flux'],xerr=dict_one[0]['x_err'], yerr=dict_one[0]['flux_err'], label=os.path.basename(file_one), fmt='r.')

        # if not (file_one.endswith(".lc")):
        #     # f1_ax1.scatter(data_array_one['x'], data_array_one['sensitivity'], marker="o",s=7,color="r" )
        #     #
        #     f1_ax1.plot(data_array_one['x'], data_array_one['sensitivity'], color='r', marker='o',  linewidth=1,linestyle='dashed', markersize=3,label=os.path.basename(file_one)+" (4 sigma sensitivity)")
        #
        #     x1 = np.subtract(data_array_one['x'],data_array_one['xerr'])
        #     x2 = np.sum([data_array_one['x'],data_array_one['xerr']],axis=0)
        #     y1 = data_array_one['sensitivity']
        #     y2 = data_array_one['sensitivity']
        #     f1_ax1.plot([x1,x2],[y1,y2],color = 'r',linestyle="dashed",  linewidth=1)

        #plot fixed flux
        if(fixed_flux!=-1):
            plt.plot([np.amin(data_array_one['x'])-np.amax(data_array_one['xerr']),np.amax(data_array_one['x'])+np.amax(data_array_one['xerr'])], [fixed_flux,fixed_flux], color='g', linestyle='dashed' , linewidth=1.5,label="Fixed Flux")

        plt.xlabel('Time MJD')
        plt.ylabel('Flux ph/cm2 s *10^-8')
        plt.ylim(0)
        #ax[0].xlabel("TT time *10^8")
        #ax[0].ylabel("Flux ph/cm2 s *10^-8")
        #ax[0].ylim(ymin=0)
        plt.grid(True)
        plt.legend(prop={'size':6})




        fig1 = plt.figure(constrained_layout=False)

        gs = GridSpec(6,1)
        f1_ax1 = fig1.add_subplot(gs[:-3,0])
        f1_ax2 = fig1.add_subplot(gs[3,0])
        f1_ax3 = fig1.add_subplot(gs[4,0])
        f1_ax4 = fig1.add_subplot(gs[5,0])

        #plt.errorbar(data_array_one['x'], data_array_one['flux'],xerr=data_array_one['xerr'], yerr=data_array_one['flux_err'], label='flux 1', fmt='r.')

        count = 0
        for detection in dict_one:
            if(count==0):
                count=1
                continue

            if(detection['sqrtts']<3):
                #add arrow
                #plt.arrow(detection['x'],detection['flux'],0,-70,width=1,head_width=10,head_starts_at_zero=True)
                f1_ax1.errorbar(detection['x'], detection['flux_ul'],xerr=detection['x_err'], yerr=0, fmt='rv')
            else:
                f1_ax1.errorbar(detection['x'], detection['flux'],xerr=detection['x_err'], yerr=detection['flux_err'], fmt='r.')

        if(dict_one[0]['sqrtts']<3):
            f1_ax1.errorbar(dict_one[0]['x'], dict_one[0]['flux_ul'],xerr=dict_one[0]['x_err'], yerr=0, label=os.path.basename(file_one), fmt='rv')
        else:
            f1_ax1.errorbar(dict_one[0]['x'], dict_one[0]['flux'],xerr=dict_one[0]['x_err'], yerr=dict_one[0]['flux_err'], label=os.path.basename(file_one), fmt='r.')

        if not (file_one.endswith(".lc")):
            # f1_ax1.scatter(data_array_one['x'], data_array_one['sensitivity'], marker="o",s=7,color="r" )
            #
            f1_ax1.plot(data_array_one['x'], data_array_one['sensitivity'], color='r', marker='o',  linewidth=1,linestyle='dashed', markersize=3,label=os.path.basename(file_one)+" (4 sigma sensitivity)")

            x1 = np.subtract(data_array_one['x'],data_array_one['xerr'])
            x2 = np.sum([data_array_one['x'],data_array_one['xerr']],axis=0)
            y1 = data_array_one['sensitivity']
            y2 = data_array_one['sensitivity']
            f1_ax1.plot([x1,x2],[y1,y2],color = 'r',linestyle="dashed",  linewidth=1)

        #plot fixed flux
        if(fixed_flux!=-1):
            f1_ax1.plot([np.amin(data_array_one['x'])-np.amax(data_array_one['xerr']),np.amax(data_array_one['x'])+np.amax(data_array_one['xerr'])], [fixed_flux,fixed_flux], color='g', linestyle='dashed' , linewidth=1.5,label="Fixed Flux")




        #LC sqrts
        f1_ax2.errorbar(data_array_one['x'], data_array_one['sqrtts'], xerr=data_array_one['xerr'],label=os.path.basename(file_one), fmt='r.')
        f1_ax2.set(xlabel="Time MJD", ylabel="sqrt(TS)")
        f1_ax2.grid(True)

        sqrtts_max = np.max(data_array_one['sqrtts'])

        if(sqrtts_max<10):
            f1_ax2.yaxis.set_major_locator(MultipleLocator(1))
        else:
            f1_ax2.yaxis.set_major_locator(MultipleLocator(round(sqrtts_max/10)))

        #LC EXP
        f1_ax3.errorbar(data_array_one['x'], data_array_one['exp'],xerr=data_array_one['xerr'], linestyle="-" , label=os.path.basename(file_one), fmt='r.')
        f1_ax3.set(xlabel="Time MJD", ylabel="Exp 10^6")
        f1_ax3.grid(True)

        #LC COUNT
        f1_ax4.errorbar(data_array_one['x'], data_array_one['count'],xerr=data_array_one['xerr'], yerr=data_array_one['count_err'], label=os.path.basename(file_one), fmt='r.')
        f1_ax4.set(xlabel="Time MJD", ylabel="Counts")
        f1_ax4.grid(True)

        f1_ax4.scatter(data_array_one['x'], data_array_one['count_bkg'], marker="o",s=7,color="r" )

        x1 = np.subtract(data_array_one['x'],data_array_one['xerr'])
        x2 = np.sum([data_array_one['x'],data_array_one['xerr']],axis=0)
        y1 = data_array_one['count_bkg']
        y2 = data_array_one['count_bkg']
        f1_ax4.plot([x1,x2],[y1,y2],color = 'r',linestyle="dashed", linewidth=1)

        max_one = np.max(data_array_one['count'])
        if(max_one<10):
            f1_ax4.yaxis.set_major_locator(MultipleLocator(1))
        else:
            f1_ax4.yaxis.set_major_locator(MultipleLocator(round(max_one/10)))





        f1_ax1.set(xlabel='Time MJD', ylabel='Flux ph/cm2 s *10^-8')
        f1_ax1.set_ylim(ymin=0)
        #ax[0].xlabel("TT time *10^8")
        #ax[0].ylabel("Flux ph/cm2 s *10^-8")
        #ax[0].ylim(ymin=0)
        f1_ax1.grid(True)
        f1_ax1.legend(prop={'size':6})


        fig1.subplots_adjust(left=0.06,bottom=0.06,right=0.97,top=0.97, wspace=0, hspace=0)
        #fig1.subplots_adjust( wspace=0, hspace=1)

        fig2, ax = plt.subplots(2,5)
        fig2.suptitle('File 1 (red): '+os.path.basename(file_one))
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
        #ax[1,3].set_title("Exp")
        ax[1,3].set(xlabel='Exp 10^6')

        if(file_one.endswith(".ap3")):
            ax[0,4].scatter(data_array_one['exp_norm'],data_array_one['flux'],color='red',s=10)
            ax[0,4].set(xlabel='Exp Norm 10^6', ylabel='Flux ph/cm2 s *10^-8')
            ax[0,4].set_title("Exp norm vs Flux")

            n, bins, patches = ax[1,4].hist(data_array_one['exp_norm'], binsize, density=True, facecolor='r', alpha=0.6,label='exp 1')
            ax[1,4].plot(bins, norm.pdf(bins, np.mean(data_array_one['exp_norm']), np.std(data_array_one['exp_norm'])), color="black", linestyle="--", alpha=0.9)

            #ax[1,4].set_title("Exp Norm")
            ax[1,4].set(xlabel='Exp Norm 10^6')


        fig2.subplots_adjust(wspace=0.4, hspace=0.2)

        fig3, ax = plt.subplots(3,3)
        fig3.suptitle('File 1 (red): '+os.path.basename(file_one))

        bins = 40

        #LC COUNT
        ax[0,0].errorbar(data_array_one['x'], data_array_one['count'],xerr=data_array_one['xerr'], yerr=data_array_one['count_err'], label=os.path.basename(file_one), fmt='r.')
        ax[0,0].set(xlabel="Time MJD", ylabel="Counts")
        ax[0,0].grid(True)

        #LC COUNT BKG
        ax[0,1].errorbar(data_array_one['x'], data_array_one['count_bkg'],xerr=data_array_one['xerr'], yerr=data_array_one['count_bkg_err'], label=os.path.basename(file_one), fmt='r.')
        ax[0,1].set(xlabel="Time MJD", ylabel="Counts Bkg")
        ax[0,1].grid(True)

        #LC EXP
        ax[1,0].errorbar(data_array_one['x'], data_array_one['exp'],xerr=data_array_one['xerr'],  label=os.path.basename(file_one), fmt='r.')
        ax[1,0].set(xlabel="Time MJD", ylabel="Exp 10^6")
        ax[1,0].grid(True)

        #LC EXP NORM
        ax[1,1].errorbar(data_array_one['x'], data_array_one['exp_norm'],xerr=data_array_one['xerr'] , label=os.path.basename(file_one), fmt='r.')
        ax[1,1].set(xlabel="Time MJD", ylabel="Exp Norm 10^6")
        ax[1,1].grid(True)


        #LC RATE
        ax[0,2].errorbar(data_array_one['x'], data_array_one['rate'],xerr=data_array_one['xerr'],yerr=data_array_one['rate_err'],  label=os.path.basename(file_one), fmt='r.')
        ax[0,2].set(xlabel="Time MJD", ylabel="Rate")
        ax[0,2].grid(True)

        #ISTOGRAMMA COUNT

        mean = np.mean(data_array_one['count'])
        std = np.std(data_array_one['count'])
        n, bins, patches = ax[2,2].hist(data_array_one['count'], binsize, density=True, facecolor='r', alpha=0.6,label="mean: "+str(round(mean,2))+",std= "+str(round(std,2)))
        ax[2,2].plot(bins, norm.pdf(bins, mean, std), color="black", linestyle="--", alpha=0.9)
        #ax[2,2].set_title("Count")
        ax[2,2].set(xlabel='Count')
        ax[2,2].legend()



        #ISTOGRAMMA COUNT BKG
        mean = np.mean(data_array_one['count_bkg'])
        std = np.std(data_array_one['count_bkg'])
        n, bins, patches = ax[1,2].hist(data_array_one['count_bkg'], binsize, density=True, facecolor='r', alpha=0.6,label="mean: "+str(round(mean,2))+",std= "+str(round(std,2)))
        ax[1,2].plot(bins, norm.pdf(bins, mean, std), color="black", linestyle="--", alpha=0.9)
        #ax[1,2].set_title("Count")
        ax[1,2].set(xlabel='Count')
        ax[1,2].legend()

        fig3.subplots_adjust(wspace=0.4, hspace=0.3)


        #plt.tight_layout()

        plt.show()

    if(mode=="2"):

        dict_two = extract_data(file_two,tstart_window_mjd,tstop_window_mjd)

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


        #plt.figure(1)
        #fig, ax = plt.subplots(ncols=1, nrows=4, constrained_layout=True)
        fig0 = plt.figure(constrained_layout=False)

        count = 0
        for detection in dict_one:
            if(count==0):
                count=1
                continue

            if(detection['sqrtts']<3):
                #add arrow
                #plt.arrow(detection['x'],detection['flux'],0,-70,width=1,head_width=10,head_starts_at_zero=True)
                plt.errorbar(detection['x'], detection['flux_ul'],xerr=detection['x_err'], yerr=0, fmt='rv')
            else:
                plt.errorbar(detection['x'], detection['flux'],xerr=detection['x_err'], yerr=detection['flux_err'], fmt='r.')

        count = 0
        for detection in dict_two:
            if(count==0):
                count=1
                continue

            if(detection['sqrtts']<3):
                #add arrow
                #plt.arrow(detection['x'],detection['flux'],0,-70,width=1,head_width=10,head_starts_at_zero=True)
                plt.errorbar(detection['x'], detection['flux_ul'],xerr=detection['x_err'], yerr=0, fmt='bv')
            else:
                plt.errorbar(detection['x'], detection['flux'],xerr=detection['x_err'], yerr=detection['flux_err'], fmt='b.')


        #plt.errorbar(data_array_one['x'], data_array_one['flux'],xerr=data_array_one['xerr'], yerr=data_array_one['flux_err'], label=os.path.basename(file_one), fmt='r.')

        if(dict_one[0]['sqrtts']<3):
            plt.errorbar(dict_one[0]['x'], dict_one[0]['flux_ul'],xerr=dict_one[0]['x_err'], yerr=0, label=os.path.basename(file_one), fmt='rv')
        else:
            plt.errorbar(dict_one[0]['x'], dict_one[0]['flux'],xerr=dict_one[0]['x_err'], yerr=dict_one[0]['flux_err'], label=os.path.basename(file_one), fmt='r.')


        if(dict_two[0]['sqrtts']<3):
            plt.errorbar(dict_two[0]['x'], dict_two[0]['flux_ul'],xerr=dict_two[0]['x_err'], yerr=0, label=os.path.basename(file_two), fmt='bv')
        else:
            plt.errorbar(dict_two[0]['x'], dict_two[0]['flux'],xerr=dict_two[0]['x_err'], yerr=dict_two[0]['flux_err'], label=os.path.basename(file_two), fmt='b.')

        # if not (file_one.endswith(".lc")):
        #     # f1_ax1.scatter(data_array_one['x'], data_array_one['sensitivity'], marker="o",s=7,color="r" )
        #     #
        #
        #     #,xerr=dict_one[0]['x_err'], yerr=0
        #     plt.errorbar(data_array_one['x'], data_array_one['sensitivity'] , color='r', marker='o',  linestyle='dashed', linewidth=1, markersize=3,label=os.path.basename(file_one)+" (4 sigma sensitivity)")
        #     x1 = np.subtract(data_array_one['x'],data_array_one['xerr'])
        #     x2 = np.sum([data_array_one['x'],data_array_one['xerr']],axis=0)
        #     y1 = data_array_one['sensitivity']
        #     y2 = data_array_one['sensitivity']
        #     plt.plot([x1,x2],[y1,y2],color = 'r',linestyle="dashed", linewidth=1)
        #
        # if not (file_two.endswith(".lc")):
        #     # f1_ax1.scatter(data_array_two['x'], data_array_two['sensitivity'], marker="o",s=7,color="b" )
        #     #
        #
        #     #,xerr=dict_two[0]['x_err'], yerr=0,
        #     plt.errorbar(data_array_two['x'], data_array_two['sensitivity'],color='b', marker='o',  linestyle='dashed', linewidth=1, markersize=3,label=os.path.basename(file_two)+" (4 sigma sensitivity)")
        #     x1 = np.subtract(data_array_two['x'],data_array_two['xerr'])
        #     x2 = np.sum([data_array_two['x'],data_array_two['xerr']],axis=0)
        #     y1 = data_array_two['sensitivity']
        #     y2 = data_array_two['sensitivity']
        #     plt.plot([x1,x2],[y1,y2],color = 'b',linestyle="dashed", linewidth=1)
        # #plot fixed flux

        if(fixed_flux!=-1):
            plt.plot([np.amin(data_array_one['x'])-np.amax(data_array_one['xerr']),np.amax(data_array_one['x'])+np.amax(data_array_one['xerr'])], [fixed_flux,fixed_flux], color='g',linestyle='dashed', linewidth=1.5,label="Fixed Flux")


        plt.xlabel('Time MJD')
        plt.ylabel('Flux ph/cm2 s *10^-8')
        plt.ylim(0)
        #ax[0].xlabel("TT time *10^8")
        #ax[0].ylabel("Flux ph/cm2 s *10^-8")
        #ax[0].ylim(ymin=0)
        plt.grid(True)
        plt.legend(prop={'size':6})


        fig1 = plt.figure(constrained_layout=False)

        gs = GridSpec(6,1)
        f1_ax1 = fig1.add_subplot(gs[:-3,0])
        #f3_ax1.set_title('gs[0, :]')
        f1_ax2 = fig1.add_subplot(gs[3,0])
        #f3_ax2.set_title('gs[0, :]')
        f1_ax3 = fig1.add_subplot(gs[4,0])
        #f3_ax3.set_title('gs[0, :]')
        f1_ax4 = fig1.add_subplot(gs[5,0])
        #f3_ax4.set_title('gs[0, :]')

        count = 0
        for detection in dict_one:
            if(count==0):
                count=1
                continue

            if(detection['sqrtts']<3):
                #add arrow
                #plt.arrow(detection['x'],detection['flux'],0,-70,width=1,head_width=10,head_starts_at_zero=True)
                f1_ax1.errorbar(detection['x'], detection['flux_ul'],xerr=detection['x_err'], yerr=0, fmt='rv')
            else:
                f1_ax1.errorbar(detection['x'], detection['flux'],xerr=detection['x_err'], yerr=detection['flux_err'], fmt='r.')

        count = 0
        for detection in dict_two:
            if(count==0):
                count=1
                continue

            if(detection['sqrtts']<3):
                #add arrow
                #plt.arrow(detection['x'],detection['flux'],0,-70,width=1,head_width=10,head_starts_at_zero=True)
                f1_ax1.errorbar(detection['x'], detection['flux_ul'],xerr=detection['x_err'], yerr=0, fmt='bv')
            else:
                f1_ax1.errorbar(detection['x'], detection['flux'],xerr=detection['x_err'], yerr=detection['flux_err'], fmt='b.')


        #plt.errorbar(data_array_one['x'], data_array_one['flux'],xerr=data_array_one['xerr'], yerr=data_array_one['flux_err'], label=os.path.basename(file_one), fmt='r.')

        if(dict_one[0]['sqrtts']<3):
            f1_ax1.errorbar(dict_one[0]['x'], dict_one[0]['flux_ul'],xerr=dict_one[0]['x_err'], yerr=0, label=os.path.basename(file_one), fmt='rv')
        else:
            f1_ax1.errorbar(dict_one[0]['x'], dict_one[0]['flux'],xerr=dict_one[0]['x_err'], yerr=dict_one[0]['flux_err'], label=os.path.basename(file_one), fmt='r.')


        if(dict_two[0]['sqrtts']<3):
            f1_ax1.errorbar(dict_two[0]['x'], dict_two[0]['flux_ul'],xerr=dict_two[0]['x_err'], yerr=0, label=os.path.basename(file_two), fmt='bv')
        else:
            f1_ax1.errorbar(dict_two[0]['x'], dict_two[0]['flux'],xerr=dict_two[0]['x_err'], yerr=dict_two[0]['flux_err'], label=os.path.basename(file_two), fmt='b.')

        if not (file_one.endswith(".lc")):
            # f1_ax1.scatter(data_array_one['x'], data_array_one['sensitivity'], marker="o",s=7,color="r" )
            #

            #,xerr=dict_one[0]['x_err'], yerr=0
            f1_ax1.errorbar(data_array_one['x'], data_array_one['sensitivity'] , color='r', marker='o',  linestyle='dashed', linewidth=1, markersize=3,label=os.path.basename(file_one)+" (4 sigma sensitivity)")
            x1 = np.subtract(data_array_one['x'],data_array_one['xerr'])
            x2 = np.sum([data_array_one['x'],data_array_one['xerr']],axis=0)
            y1 = data_array_one['sensitivity']
            y2 = data_array_one['sensitivity']
            f1_ax1.plot([x1,x2],[y1,y2],color = 'r',linestyle="dashed", linewidth=1)

        if not (file_two.endswith(".lc")):
            # f1_ax1.scatter(data_array_two['x'], data_array_two['sensitivity'], marker="o",s=7,color="b" )
            #

            #,xerr=dict_two[0]['x_err'], yerr=0,
            f1_ax1.errorbar(data_array_two['x'], data_array_two['sensitivity'],color='b', marker='o',  linestyle='dashed', linewidth=1, markersize=3,label=os.path.basename(file_two)+" (4 sigma sensitivity)")
            x1 = np.subtract(data_array_two['x'],data_array_two['xerr'])
            x2 = np.sum([data_array_two['x'],data_array_two['xerr']],axis=0)
            y1 = data_array_two['sensitivity']
            y2 = data_array_two['sensitivity']
            f1_ax1.plot([x1,x2],[y1,y2],color = 'b',linestyle="dashed", linewidth=1)
        #plot fixed flux

        if(fixed_flux!=-1):
            f1_ax1.plot([np.amin(data_array_one['x'])-np.amax(data_array_one['xerr']),np.amax(data_array_one['x'])+np.amax(data_array_one['xerr'])], [fixed_flux,fixed_flux], color='g',linestyle='dashed', linewidth=1.5,label="Fixed Flux")


        f1_ax1.set(xlabel='Time MJD', ylabel='Flux ph/cm2 s *10^-8')
        f1_ax1.set_ylim(ymin=0)
        #ax[0].xlabel("TT time *10^8")
        #ax[0].ylabel("Flux ph/cm2 s *10^-8")
        #ax[0].ylim(ymin=0)
        f1_ax1.grid(True)
        f1_ax1.legend(prop={'size':6})

        #LC TS
        f1_ax2.errorbar(data_array_one['x'], data_array_one['sqrtts'],xerr=data_array_one['xerr'], label=os.path.basename(file_one), fmt='r.')
        f1_ax2.errorbar(data_array_two['x'], data_array_two['sqrtts'],xerr=data_array_two['xerr'], label=os.path.basename(file_two), fmt='b.')
        f1_ax2.set(xlabel="Time MJD", ylabel="sqrt(TS)")
        f1_ax2.grid(True)

        min_one = np.amin(data_array_one['sqrtts'])
        min_two = np.amin(data_array_two['sqrtts'])

        max_one = np.max(data_array_one['sqrtts'])
        max_two = np.max(data_array_two['sqrtts'])

        sqrtts_min = np.amin([min_one,min_two])
        sqrtts_max = np.amax([max_one,max_two])

        #grid_array = np.arange(int(sqrtts_min),int(sqrtts_max),1)

        #for grid_line in grid_array.tolist():
        #    f1_ax2.axhline(grid_line, color='gray', linewidth=0.5)
        if(sqrtts_max<10):
            f1_ax2.yaxis.set_major_locator(MultipleLocator(1))
        else:
            f1_ax2.yaxis.set_major_locator(MultipleLocator(round(sqrtts_max/10)))
        #f1_ax2.set_axisbelow(True)

        #LC EXP
        f1_ax3.errorbar(data_array_one['x'], data_array_one['exp'],xerr=data_array_one['xerr'], linestyle="-" ,  label=os.path.basename(file_one), fmt='r.')
        f1_ax3.errorbar(data_array_two['x'], data_array_two['exp'],xerr=data_array_two['xerr'], linestyle="-" ,  label=os.path.basename(file_two), fmt='b.')
        f1_ax3.set(xlabel="Time MJD", ylabel="Exp 10^6")
        f1_ax3.grid(True)

        #LC COUNT
        f1_ax4.errorbar(data_array_one['x'], data_array_one['count'],xerr=data_array_one['xerr'], yerr=data_array_one['count_err'], label=os.path.basename(file_one), fmt='r.')
        f1_ax4.errorbar(data_array_two['x'], data_array_two['count'],xerr=data_array_two['xerr'], yerr=data_array_two['count_err'], label=os.path.basename(file_two), fmt='b.')

        f1_ax4.scatter(data_array_one['x'], data_array_one['count_bkg'], marker="o",s=7,color="r" )
        f1_ax4.scatter(data_array_two['x'], data_array_two['count_bkg'], marker="o",s=7,color="b" )

        x1 = np.subtract(data_array_one['x'],data_array_one['xerr'])
        x2 = np.sum([data_array_one['x'],data_array_one['xerr']],axis=0)
        y1 = data_array_one['count_bkg']
        y2 = data_array_one['count_bkg']
        f1_ax4.plot([x1,x2],[y1,y2],color = 'r',linestyle="dashed",linewidth=1,)

        x1 = np.subtract(data_array_two['x'],data_array_two['xerr'])
        x2 = np.sum([data_array_two['x'],data_array_two['xerr']],axis=0)
        y1 = data_array_two['count_bkg']
        y2 = data_array_two['count_bkg']
        f1_ax4.plot([x1,x2],[y1,y2],color = 'b',linestyle="dashed",linewidth=1,)

        f1_ax4.set(xlabel="Time MJD", ylabel="Counts")
        f1_ax4.grid(True)

        min_one = np.amin(data_array_one['count'])
        min_two = np.amin(data_array_two['count'])

        max_one = np.max(data_array_one['count'])
        max_two = np.max(data_array_two['count'])

        cts_min = np.min([min_one,min_two])
        cts_max = np.max([max_one,max_two])

        #grid_array = np.arange(int(cts_min),int(cts_max),1)

        #for grid_line in grid_array.tolist():
        #    f1_ax4.axhline(grid_line, color='gray', linewidth=0.5)
        if(cts_max<10):
            f1_ax4.yaxis.set_major_locator(MultipleLocator(1))
        else:
            f1_ax4.yaxis.set_major_locator(MultipleLocator(round(cts_max/10)))

        #fig1.tight_layout()
        fig1.subplots_adjust(left=0.06,bottom=0.06,right=0.97,top=0.97,wspace=0,hspace=0)



        fig2, ax = plt.subplots(2,5)
        fig2.suptitle('File 1 (red): '+os.path.basename(file_one)+', File 2 (blue): '+os.path.basename(file_two))
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
            #ax[1,4].set_title("Exp Norm")
            ax[1,4].set(xlabel='Exp Norm 10^6')

        if(file_two.endswith(".ap3")):

            ax[0,4].scatter(data_array_two['exp_norm'],data_array_two['flux'],color='blue',s=10)
            ax[0,4].set(xlabel='Exp Norm 10^6', ylabel='Flux ph/cm2 s *10^-8')
            ax[0,4].set_title("Exp Norm vs Flux")

            n, bins, patches = ax[1,4].hist(data_array_two['exp_norm'], binsize, density=True, facecolor='b', alpha=0.6,label='exp 1')
            ax[1,4].plot(bins, norm.pdf(bins, np.mean(data_array_two['exp_norm']), np.std(data_array_two['exp_norm'])), color="black", linestyle="--", alpha=0.9)
            #ax[1,4].set_title("Exp Norm")
            ax[1,4].set(xlabel='Exp Norm 10^6')

        mean = np.mean(data_array_one['exp'])
        std = np.std(data_array_one['exp'])
        n, bins, patches = ax[1,3].hist(data_array_one['exp'], binsize, density=True, facecolor='r', alpha=0.6,label="mean: "+str(round(mean,2))+",std= "+str(round(std,2)))
        ax[1,3].plot(bins, norm.pdf(bins, mean, std), color="black", linestyle="--", alpha=0.9)

        mean = np.mean(data_array_two['exp'])
        std = np.std(data_array_two['exp'])
        n, bins, patches = ax[1,3].hist(data_array_two['exp'], binsize, density=True, facecolor='b', alpha=0.6,label="mean: "+str(round(mean,2))+",std= "+str(round(std,2)))
        ax[1,3].plot(bins, norm.pdf(bins, mean, std), color="black", linestyle="--", alpha=0.9)
        #ax[1,3].set_title("Exp")
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

        fig2.subplots_adjust(wspace=0.4, hspace=0.2)



        fig3, ax = plt.subplots(3,3)
        fig3.suptitle('File 1 (red): '+os.path.basename(file_one)+', File 2 (blue): '+os.path.basename(file_two))

        bins = 40

        #LC COUNT
        ax[0,0].errorbar(data_array_one['x'], data_array_one['count'],xerr=data_array_one['xerr'], yerr=data_array_one['count_err'], label=os.path.basename(file_one), fmt='r.')
        ax[0,0].errorbar(data_array_two['x'], data_array_two['count'],xerr=data_array_two['xerr'], yerr=data_array_two['count_err'], label=os.path.basename(file_two), fmt='b.')
        ax[0,0].set(xlabel="Time MJD", ylabel="Counts")
        ax[0,0].grid(True)

        #LC COUNT BKG
        ax[0,1].errorbar(data_array_one['x'], data_array_one['count_bkg'],xerr=data_array_one['xerr'], yerr=data_array_one['count_bkg_err'], label=os.path.basename(file_one), fmt='r.')
        ax[0,1].errorbar(data_array_two['x'], data_array_two['count_bkg'],xerr=data_array_two['xerr'], yerr=data_array_two['count_bkg_err'], label=os.path.basename(file_two), fmt='b.')
        ax[0,1].set(xlabel="Time MJD", ylabel="Counts Bkg")
        ax[0,1].grid(True)

        #LC EXP
        ax[1,0].errorbar(data_array_one['x'], data_array_one['exp'],xerr=data_array_one['xerr'],  label=os.path.basename(file_one), fmt='r.')
        ax[1,0].errorbar(data_array_two['x'], data_array_two['exp'],xerr=data_array_two['xerr'],  label=os.path.basename(file_two), fmt='b.')
        ax[1,0].set(xlabel="Time MJD", ylabel="Exp 10^6")
        ax[1,0].grid(True)

        #LC EXP NORM
        ax[1,1].errorbar(data_array_one['x'], data_array_one['exp_norm'],xerr=data_array_one['xerr'] , label=os.path.basename(file_one), fmt='r.')
        ax[1,1].errorbar(data_array_two['x'], data_array_two['exp_norm'],xerr=data_array_two['xerr'], label=os.path.basename(file_two), fmt='b.')
        ax[1,1].set(xlabel="Time MJD", ylabel="Exp Norm 10^6")
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
        ax[0,2].set(xlabel="Time MJD", ylabel="Rate")
        ax[0,2].grid(True)

        #ISTOGRAMMA COUNT

        mean = np.mean(data_array_one['count'])
        std = np.std(data_array_one['count'])
        n, bins, patches = ax[2,2].hist(data_array_one['count'], binsize, density=True, facecolor='r', alpha=0.6,label="mean: "+str(round(mean,2))+",std= "+str(round(std,2)))
        ax[2,2].plot(bins, norm.pdf(bins, mean, std), color="black", linestyle="--", alpha=0.9)
        #ax[2,2].set_title("Count")
        ax[2,2].set(xlabel='Count')
        ax[2,2].legend()


        mean = np.mean(data_array_two['count'])
        std = np.std(data_array_two['count'])
        n, bins, patches = ax[2,2].hist(data_array_two['count'], binsize, density=True, facecolor='b', alpha=0.6,label="mean: "+str(round(mean,2))+",std= "+str(round(std,2)))
        ax[2,2].plot(bins, norm.pdf(bins, mean, std), color="black", linestyle="--", alpha=0.9)
        #ax[2,2].set_title("Count")
        ax[2,2].set(xlabel='Count')
        ax[2,2].legend()


        #ISTOGRAMMA COUNT BKG
        mean = np.mean(data_array_one['count_bkg'])
        std = np.std(data_array_one['count_bkg'])
        n, bins, patches = ax[1,2].hist(data_array_one['count_bkg'], binsize, density=True, facecolor='r', alpha=0.6,label="mean: "+str(round(mean,2))+",std= "+str(round(std,2)))
        ax[1,2].plot(bins, norm.pdf(bins, mean, std), color="black", linestyle="--", alpha=0.9)
        #ax[1,2].set_title("Count")
        ax[1,2].set(xlabel='Count')
        ax[1,2].legend()


        mean = np.mean(data_array_two['count_bkg'])
        std = np.std(data_array_two['count_bkg'])
        n, bins, patches = ax[1,2].hist(data_array_two['count_bkg'], binsize, density=True, facecolor='b', alpha=0.6,label="mean: "+str(round(mean,2))+",std= "+str(round(std,2)))
        ax[1,2].plot(bins, norm.pdf(bins, mean, std), color="black", linestyle="--", alpha=0.9)
        #ax[1,2].set_title("Count")
        ax[1,2].set(xlabel='Count')
        ax[1,2].legend()

        fig3.subplots_adjust(wspace=0.4, hspace=0.3)

        plt.show()

if __name__ == "__main__":
    main()