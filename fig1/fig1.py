# uses downloaded sac files and inventory (fiji_sac_get.py) to create Figure 1 of Taup 3.2, 2026, paper.
# comment the sys.exit() around line 143 to plot 1B.
# Shubh Agrawal
# USC, SC, July 2026
from obspy import read, Stream, UTCDateTime,read_events
from obspy.core.event import Origin, Catalog
from obspy.core.inventory.inventory import read_inventory
import numpy as np
from obspy.clients.fdsn import Client
from obspy.clients.fdsn import RoutingClient
import datetime
import os
import sys
import taup
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

from fiji_sac_get import getData, inv_fiji


def get_unique_phases(jsonfinds_small):
    seen=set()
    arrivals_unique = []
    for arr in jsonfinds_small.arrivals:
        rp = arr.rayparam
        rp = round(arr.rayparam, 4)
        if rp not in seen:
            arrivals_unique.append(arr)
            seen.add(rp)
    findsphases=[]
    for arr in arrivals_unique:
        findsphases.append(arr.phase)
    return findsphases

def makePlots(fullImage, taupserver):

    ####
    mpl.rcParams.update({'font.size': 16.5})
    #Eq deets
    starttime= UTCDateTime('2018-08-19T00:19:40') #Fiji
    eq_lat,eq_long,eq_depth=-18.1125,-178.1530,600 #Fiji
    endtime= starttime+7000
    ###
    inventory_big=read_inventory(inv_fiji)
    # # inventory_big.plot(label=False,color_per_network=True,resolution='i',continent_fill_color='honeydew',alpha=.5)
    ### read eq data
    stream_all=read('sac_fiji_18/*.sac')
    stream_all.taper(type='cosine',max_percentage=.1)
    ### TAUP server
    modelname = 'iasp91'
    #Distance query to get garc from earthquake to station
    params = taup.DistazQuery()
    params.model(modelname)
    params.event(eq_lat,eq_long)
    params.sourcedepth(eq_depth)
    for tr in stream_all:
        inv = inventory_big.select(network=tr.stats.network,station=tr.stats.station)
        sta = inv[0][0]
        params.station(sta.latitude,sta.longitude)
        jsondists = params.calc(taupserver)
        tr.stats.distance=jsondists.distances[0].deg
    #### Find query to get phases using earthquake depth and --max
    params_f = taup.FindQuery()
    params_f.model(modelname)
    params_f.max(3)
    params_f.sourcedepth(eq_depth)
    params_f.exclude('20,210,moho,410,660')
    find_max3 = params_f.calc(taupserver)
    find_max3.foundphases = list(dict.fromkeys(find_max3.foundphases))
    print('length of phases in entire window:',len(find_max3.foundphases))

    ### taup find query for time between 2000,3500 and dist 80,125
    # cmdLine = params.asCommandLine(taupserver)
    # taupserver.verbose = True
    params_f.max(5)
    params_f.time([1900,2900])
    params_f.deg([98])
    params_f.rayparamdeg([7,25])
    jsonfinds_small = params_f.calc(taupserver)
    cmdLine = params_f.asCommandLine(taupserver)
    print('command line prompt:',cmdLine)

    findsphases_subset=get_unique_phases(jsonfinds_small)
    # print('length of found phases in smaller window:',len(jsonfinds_small.arrivals))
    print('length of unique phases in smaller time window:',len(findsphases_subset))
    print(findsphases_subset)

    # Curve query to plot traveltime-dist curve for the phases in find
    params_c = taup.CurveQuery()
    params_c.model(modelname)
    params_c.sourcedepth(eq_depth)
    params_c.phase(find_max3.foundphases)
    jsoncurve = params_c.calc(taupserver)
    ### asking curves for the subset of phases
    params_c.phase(findsphases_subset)
    jsoncurve_subset = params_c.calc(taupserver)

    ##
    stream_all.sort(['distance'],reverse=True)
    ####
    ### keep the first trace in a degree bin!
    st_bin = Stream()
    last_bin = None
    for tr in stream_all:
        d_bin = int(tr.stats.distance) # degree bin hack
        if d_bin != last_bin:
            st_bin.append(tr)
            last_bin = d_bin

    plt.ion()
    fig = plt.figure(figsize=(17, 8))
    sns.set_style("whitegrid",{"axes.facecolor": ".2","grid.color": ".6", "grid.linestyle": ":"})
    # [left, bottom, width, height]
    ax = fig.add_axes([0.07, 0.1, 0.88, 0.8])
    for tr in st_bin:#stream_all[::4]
        t = tr.times(reftime=starttime)
        ax.plot(t, 1.5*tr.data/tr.data.max() + tr.stats.distance,lw=.55,c='white')

    ax.xaxis.set_minor_locator(MultipleLocator(250))
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.yaxis.set_minor_locator(MultipleLocator(5))
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_label_position('bottom')
    ax.xaxis.tick_bottom()
    ax.set_ylabel('Distance ($^\\circ$)')
    ax.set_xlabel('Time after event (s)')
    ax.set_xlim(0,6500)
    ax.set_ylim(55,162)
    plt.title(f"Fiji August 2018 Mw 8.2")
    plt.gca().invert_yaxis()
    ax.yaxis.grid(False)
    if fullImage:

        ### this bit needs to be commented while plotting B.
        for curve in jsoncurve.curves:
            # if curve.label[0]=='p': #useful for shallow eartquakes as depth phases are not well seperated!
            #     continue
            if curve.label[-1] in ['s','S']: # phases with last leg as S wave
                for seg in curve.segments:
                    ax.plot(seg.y, seg.x,lw=.35,c='lightgreen',alpha=.75,zorder=2)
            else:
                for seg in curve.segments:
                    ax.plot(seg.y, seg.x,lw=.35,c='lightpink',alpha=.75,zorder=2)

        fig.savefig('fiji_55deg_Z_max3_.pdf')
        fig.savefig('fiji_55deg_Z_max3_.png', dpi=600, pad_inches=0.4)#bbox_inches='tight',
        # plt.show()
    else:
        ### Fig1B. sub-plot for zoom in section plot for time: 1800-3200sec!
        for curve in jsoncurve_subset.curves:
            if curve.label[-1] in ['s','S']: # phases with last leg as S wave
                for seg in curve.segments:
                    ax.plot(seg.y, seg.x,lw=.85,c='lightgreen',ls='-',alpha=.95,zorder=2)
            else:
                for seg in curve.segments:
                    ax.plot(seg.y, seg.x,lw=.85,c='lightpink',ls='-',alpha=.95,zorder=2)

        ax.set_xlim(1800,3199)
        ax.set_ylim(80,125)
        ax.xaxis.set_minor_locator(MultipleLocator(100))
        ax.xaxis.set_major_locator(MultipleLocator(200))
        plt.show()
        ##### shenanigans to plot phase name as text
        i=0
        t_last=None
        for arr in  jsonfinds_small.arrivals:
            t_c=arr.time
            if t_last==None or t_c-t_last>25: #only typing one of the phases if they arrive within 25 sec of each other.
                if arr.phase[-1] in ['s','S']:
                    ax.text(arr.time,98,arr.phase,bbox={'facecolor': 'white', 'alpha': 0.85, 'pad': 1},horizontalalignment='center', color='seagreen',fontsize=13)
                    t_last=arr.time
                else:
                    ax.text(arr.time,98-i,arr.phase,bbox={'facecolor': 'white', 'alpha': 0.85, 'pad': 1},horizontalalignment='center', color='palevioletred',fontsize=13)
                    i=i-2
        plt.gca().invert_yaxis()
        plt.title("")
        fig.savefig('fiji_55deg_Z_max5_time_.pdf')
        fig.savefig('fiji_55deg_Z_max5_time_.png', dpi=600, pad_inches=0.4)
        #
        # bin/taup find --max 3 --evdepth 600 --exclude 20,210,moho,410,660 --pwaveonly --mod iasp91 --time 2400 2600 --showrayparam --deg 100

def fig1(taupserver):
    if not os.path.exists(inv_fiji):
        getData()
    makePlots(True, taupserver) # fig a
    makePlots(False, taupserver) # fig b

def main():
    taup_path="~/Research/sct_wat/TauP/build/install/TauP/bin/taup"
    taup_path="~/Code/seis/TauP/build/install/TauP/bin/taup"

    with taup.TauPServer(taup_path=taup_path) as taupserver:
        fig1(taupserver)

if __name__ == '__main__':
    main()
