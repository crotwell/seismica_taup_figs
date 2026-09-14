## using taup curve, we compare different seismic attribitues.
# creates Figure 4 of Taup 3.2, 2026, paper.
# Shubh Agrawal
# USC, July 2026
##
import numpy as np
import taup
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys
import seaborn as sns
from collections import defaultdict
from matplotlib.colors import to_rgba
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from obspy.clients.fdsn import Client

def plot_attribute(ax,jsonCurve,cl):
    for curve in jsonCurve.curves:
        for seg in curve.segments:
            ax.plot(seg.x, seg.y,c=cl,linestyle='-',linewidth=1.5,alpha=.95,label=curve.label)
            ax.set_ylabel(curve.y)

def makePlot(taupserver):
    mpl.rcParams.update({'font.size': 16})
    ###
    client = Client("USGS")
    eq_id='us6000sasz' # Malaysia deep eq

    catalog = client.get_events(eventid=eq_id)
    origin=catalog[0].origins[0]
    eq_lat,eq_long,eq_depth=origin.latitude,origin.longitude,origin.depth/1000
    eq_mw=catalog[0].magnitudes[0].mag
    nodal_plane=catalog[0].focal_mechanisms[0].nodal_planes['nodal_plane_1']
    strike,dip,rake=nodal_plane.strike,nodal_plane.dip,nodal_plane.rake

    """
    #radian, radian180, degree, degree180, kilometer, kilometer180, rayparamrad, rayparamdeg, rayparamkm,
     time, tau, takeoffangle, incidentangle, maxdepth, dpddelta, dpddeg, amp, amppsv, ampsh, phase,
     phasepsv, phasesh, phasedeg, phasedegpsv, phasedegsh, unwrapphasedeg, unwrapphasedegpsv, unwrapphasedegsh,
     geospread, refltran, refltranpsv, refltransh, index, tstar, attenuation, energygeospread, pathlength,
     radiation, radiationpsv, radiationsh, intcaustic
    """

    yAxis = ['time',"rayparamdeg", "maxdepth",'pathlength', "amppsv", "refltran","attenuation",'radiationpsv']
    cl=['slateblue']

    vel_mod='ak135fcont'
    phase_list=['Smp','ScSP','PKKP']
    phase_list=['PP']

    plt.ion()
    fig, axs = plt.subplots(2, 4,figsize=(15, 10),sharex=False,sharey=False,constrained_layout=True)
    sns.set_style("whitegrid",{"axes.facecolor": "ghostwhite","grid.color": ".4", "grid.linestyle": ":"})
    ax1 = axs[0, 0]
    ax2 = axs[0, 1]
    ax3 = axs[0, 2]
    ax4 = axs[0, 3]
    ax5 = axs[1, 0]
    ax6 = axs[1, 1]
    ax7 = axs[1, 2]
    ax8 = axs[1, 3]
    axx=[ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8]

    params = taup.CurveQuery()
    params.model(vel_mod)
    params.mw(eq_mw)
    params.sourcedepth(eq_depth)
    params.xaxis('degree')
    params.strikediprake(strike,dip,rake)
    params.az(0)
    for j, yAtt in enumerate(yAxis):
        for i, phase in enumerate(phase_list):
            params.yaxis(yAtt)
            params.phase(phase)
            jsonCurve = params.calc(taupserver)
            cmdLine = params.asCommandLine(taupserver)
            print('command line prompt:',cmdLine)
            plot_attribute(axx[j],jsonCurve,cl[i])
    ax1.set_ylabel('Time (s)')
    ax2.set_ylabel('Rayparam (s/$^\circ$)')
    ax3.set_ylabel('Max depth (km)')
    ax4.set_ylabel('Path length (km)')
    ax5.set_ylabel('Amplitude P$_{sv}$ (m)')
    ax6.set_ylabel('Reflection-transmission coeff.')
    ax7.set_ylabel('Attenuation')
    ax8.set_ylabel('Radiation P$_{sv}$')

    fig.supxlabel('Distance ($^\circ$)', x=0.52)

    xmin, xmax = 50, 200
    for ax in axx:
        ax.set_xlim(50,200)
        for line in ax.lines:
            x = line.get_xdata()
            y = line.get_ydata()

            mask = (x >= xmin) & (x <= xmax)
            if np.any(mask):
                ax.set_ylim(np.nanmin(y[mask]), np.nanmax(y[mask]))

        ax.set_xlim(xmin, xmax)

    #ax8.legend()
    plt.savefig('attributes_july28.pdf')
    plt.savefig('attributes_july28.png',dpi=400,bbox_inches='tight', pad_inches=0.1)
    ###


def fig4(taupserver):
    makePlot(taupserver)

def main():
    taup_path="~/Research/sct_wat/TauP/build/install/TauP/bin/taup"
    taup_path="~/Code/seis/TauP/build/install/TauP/bin/taup"

    with taup.TauPServer(taup_path=taup_path) as taupserver:
        fig4(taupserver)

if __name__ == '__main__':
    main()
