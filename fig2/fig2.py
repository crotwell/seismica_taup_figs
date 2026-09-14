# using the taup's amplitude feature, compares amplitude of core phases!
# generates Figure 2 of Taup 3.2, 2026, paper. Use 'plot_A' switch to either plot A or B.
# Shubh Agrawal
# USC, July 2026

import numpy as np
import taup
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys
import seaborn as sns
from obspy import read, Stream, UTCDateTime,read_events
from obspy.geodetics import gps2dist_azimuth
from collections import defaultdict
from matplotlib.colors import to_rgba
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

###
def get_dict_amps(TimeResult):
    """
    reads taup_py output and gets phase names and amplitudes. saves highest amp in case of multiple arrival at a distance..
    """
    grouped = defaultdict()
    for a in TimeResult.arrivals:
        phase=a.phase
        amp=float(a.amp.factorpsv)
        prev = grouped.get(phase)
        # bit for triplication cases
        if prev is None or abs(amp) > abs(prev):
            if prev != None: #debug step
                pass
                # print(f"prev:{prev}; phase:{phase}; abs amp{abs(amp)}")
            grouped[phase] = abs(amp)

    return grouped

def makePlot(plot_A, taupserver):
    mpl.rcParams.update({'font.size': 14.5})
    # https://earthquake.usgs.gov/earthquakes/eventpage/us1000gcii/executive
    # https://earthquake.usgs.gov/earthquakes/eventpage/us6000m52p
    eventdepth=607
    # phases=['PKP','PKIKP','SKIKS','PKJKP','SKJKS']
    phases=['PKP','PKIKP','PKJKP']

    plt.ion()
    plt.figure(figsize=(14, 5))
    ax = plt.axes()
    ax.set_facecolor(to_rgba('darkseagreen', alpha=0.1))
    plt.style.use('ggplot')
    dists = []
    ratios = []
    amps_J=[]
    amps_I=[]
    amps_K=[]
    colors=['mediumpurple','cadetblue','indianred','skyblue','darkgrey']
    params = taup.TimeQuery()
    params.model('prem')
    params.sourcedepth(eventdepth)
    if plot_A:
        params.az(80)
        markr='*'
    else:
        dist=210
        params.degree(dist)
        markr='X'

    params.amp(True)
    params.mw(6.6)
    params.strikediprake(17,7,-62)# 17°	7° -62 Brazio
    phase_ratios = defaultdict(list)
    for i,phase in enumerate(phases):
        params.phase(phase)
        j=0
        if plot_A:
            for dist in np.arange(60.0,250,2.5):
                params.degree([dist])
                TimeResult = params.calc(taupserver)
                amps = get_dict_amps(TimeResult)
                if len(amps)== 0:
                    continue
                if amps[phase]!=0:
                    if j==0:
                        plt.scatter(dist, amps[phase], marker=markr, alpha=.8,s=45, color=colors[i],zorder=10,label=phase)
                        j=+1
                    else:
                        plt.scatter(dist, amps[phase], marker=markr, alpha=.8,s=45, color=colors[i],zorder=10)
        else:
            for az in np.arange(0.0,360,6):
                params.az([az])
                TimeResult = params.calc(taupserver)
                amps = get_dict_amps(TimeResult)
                if len(amps)== 0:
                    continue
                if amps[phase]!=0:
                    if j==0:
                        plt.scatter(az, amps[phase], marker=markr, alpha=.8,s=45, color=colors[i],zorder=10,label=phase)
                        j=+1
                    else:
                        plt.scatter(az, amps[phase], marker=markr, alpha=.8,s=45, color=colors[i],zorder=10)

    ax.set_yscale("log")
    ax.xaxis.set_minor_locator(MultipleLocator(20))
    ax.xaxis.set_major_locator(MultipleLocator(40))
    # plt.legend(loc='upper left',fontsize='14')
    plt.legend(loc='lower right',fontsize='15')

    plt.ylabel("Amplitude ($P_{sv}$)")#PKIKP/ PKJKP
    if plot_A:
        plt.xlabel("Distance ($^\\circ$)")
        ax.axvline(x=210,ls='--',lw=1.5,c='darkgrey',zorder=1)
    else:
        plt.xlabel("Azimuth ($^\\circ$)")
        ax.axvline(x=80,ls='--',lw=1.5,c='darkgrey',zorder=1)
    plt.tight_layout()

    # sys.exit()
    if plot_A:
        plt.savefig('mx6.6_3phases_az80.pdf')
        plt.savefig('mx6.6_3phases_az80.png',dpi=400,bbox_inches='tight', pad_inches=0.1)
    else:
        plt.savefig('mx6.6_3phases_dist_210.pdf')
        plt.savefig('mx6.6_3phases_dist_210.png',dpi=400,bbox_inches='tight', pad_inches=0.1)

    # plt.show()


def fig2(taupserver):
    makePlot(True, taupserver)
    makePlot(False, taupserver)


def main():
    taup_path="~/Research/sct_wat/TauP/build/install/TauP/bin/taup"
    taup_path="~/Code/seis/TauP/build/install/TauP/bin/taup"

    with taup.TauPServer(taup_path=taup_path) as taupserver:
        fig2(taupserver)

if __name__ == '__main__':
    main()
