# compares geodetic, geocentric, and spherical distances for an earthq and a grid of stations around it.
##
import numpy as np
import taup
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys
from matplotlib.colors import to_rgba
from matplotlib.lines import Line2D
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
from scipy.interpolate import griddata
import matplotlib.tri as tri
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import TwoSlopeNorm,CenteredNorm
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from ellipticipy import ellipticity_correction
from obspy.taup import TauPyModel
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees, locations2degrees

###

mpl.rcParams.update({'font.size': 15})
###
def spherical_distance(lat1, lon1, lat2, lon2):
    # lat1, lat2 = np.radians([lat1, lat2])
    lat1 = np.radians(lat1)
    lat2 = np.radians(lat2)

    dlat = lat2 - lat1
    dlon = np.radians(lon2 - lon1)

    a = np.sin(dlat/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2)**2
    return 2 * np.degrees(np.arcsin(np.sqrt(a)))

def make_station_grid(eq_lat, eq_lon, spacing=5, max_dist=50):

    lats = np.arange(-90, 90.1, spacing)
    lons = np.arange(-180, 180.1, spacing)
    lon2d, lat2d = np.meshgrid(lons, lats)

    # spherical GRC from earthquake
    dist_deg = spherical_distance(np.full(lat2d.shape, eq_lat),np.full(lon2d.shape, eq_lon),lat2d,lon2d)

    mask = (dist_deg <= max_dist) & ( spacing <=dist_deg)

    stations = np.column_stack([lat2d[mask],lon2d[mask]])
    # mask = [dist_deg <= max_dist,spacing <=dist_deg]

    return stations
def plot_stations(eq_lat, eq_lon, stations):
    fig = plt.figure(figsize=(12, 7))
    ax = plt.axes(projection=ccrs.Robinson())
    ax.set_global()
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    # ax.add_feature(cfeature.BORDERS, linewidth=0.3)

    ax.scatter(stations[:, 1],stations[:, 0],s=5,marker='^',color='teal',transform=ccrs.PlateCarree(),alpha=.65)
    ax.scatter(eq_lon,eq_lat,s=90,marker='*',color='darkred',edgecolor='darkred',transform=ccrs.PlateCarree(),zorder=5)
    plt.show()

def plot_dist_diff(station_data,param1='dist_geod',param2='dist_sph',eq_lon=0,eq_lat=45,cmap='cividis',label="$ \Delta $ Time (s)",tick_space=None,figname='map.png'):
    data = pd.DataFrame.from_dict(station_data)

    if param1=='time_correction':
        data['param_diff']=data[param1]
    else:
        data['param_diff']=data[param1]-data[param2]
    max_val = np.max(np.abs(data['param_diff']))
    print(f"max/min val for {param1} - {param2}: {np.max(data['param_diff']),np.min(data['param_diff'])}..\n")

    triang = tri.Triangulation(data['lon'],data['lat'])
    ##
    fig = plt.figure(figsize=(14, 7))
    ax = plt.axes(projection=ccrs.Robinson())
    ax.set_global()

    levels = np.linspace(-max_val, max_val, 51)

    cf = ax.tricontourf(triang,data['param_diff'],cmap=cmap,levels=levels,transform=ccrs.PlateCarree())#,extend='both')

    # ax.coastlines(linewidth=0.8)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)

    cb = plt.colorbar(cf,ax=ax,orientation='horizontal',pad=0.05,shrink=0.5,fraction=.15)
    ax.scatter(eq_lon,eq_lat,s=130,marker='*',color='royalblue',edgecolor='royalblue',transform=ccrs.PlateCarree(),zorder=5)

    # cb.ax.set_yticklabels(["{:.1}".format(i) for i in cb.get_ticks()]) # ticks
    if tick_space!=None:
        cb.ax.xaxis.set_major_locator(MultipleLocator(tick_space))
    # cb.ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    # cb.ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    cb.set_label(label)
    plt.savefig(figname,dpi=600,bbox_inches='tight', pad_inches=0.1)
    plt.show()

def spherical_dist_azi(lat_s, lon_s, lat_r, lon_r):
    sph_dist_deg = locations2degrees(lat_s, lon_s, lat_r, lon_r) #obspy haversine

    phi_s = np.radians(lat_s)
    phi_r = np.radians(lat_r)
    dlon = np.radians(lon_r - lon_s)

    azimuth = np.degrees(np.arctan2(np.sin(dlon) * np.cos(phi_r),
        np.cos(phi_s) * np.sin(phi_r)
        - np.sin(phi_s) * np.cos(phi_r) * np.cos(dlon)))

    azimuth %= 360

    return sph_dist_deg, azimuth

def ellipticipy_py_corr(phase=['SS'],spacing=5, max_dist=50,eq_lon=0,eq_lat=45):
    """
    Uses ellipticipy (https://github.com/StuartJRussell/EllipticiPy/tree/master/src)
    to calculate time correction at geodetic distances calculated using obspy..
    """
    model = TauPyModel('iasp91')
    stations = make_station_grid(eq_lat, eq_lon, spacing=spacing, max_dist=max_dist)
    eventdepth=30
    station_data_epy=[]
    for st in stations:
         dist_sp_deg,azi_sph=spherical_dist_azi(eq_lat, eq_lon, st[0], st[1])
         # dist_m, az, baz = gps2dist_azimuth(eq_lat, eq_lon, st[0], st[1])
         # dist_d=kilometers2degrees(dist_m / 1000.)
         arrivals = model.get_ray_paths(source_depth_in_km = eventdepth, distance_in_degree = dist_sp_deg, phase_list = phase)
         correction=ellipticity_correction(arrivals, azimuth = azi_sph, source_latitude = eq_lat)
         # print(f" correction val:  {correction}, for dist={dist_d}")

         station_data_epy.append({
             'lat': st[0],
             'lon': st[1],
             'dist_sp_deg': dist_sp_deg,
             'time_sph': arrivals[0].time,
             'time_correction': correction[0],
             'phase': phase})

    print(f"lenth of station data ePy: {len(station_data_epy)}")
    return station_data_epy,arrivals
###
def cal_dist_time(taupserver,phase=['P'],spacing=5, max_dist=50,eq_lon=0,eq_lat=45):
    phases=(phase)
    eventdepth=30
    # eq_lat,eq_lon=30,0
    evt=(eq_lat,eq_lon)
    params = taup.TimeQuery()
    params.model('iasp91')
    params.sourcedepth(eventdepth)
    params.event(eq_lat,eq_lon)
    stations = make_station_grid(eq_lat, eq_lon, spacing=spacing, max_dist=max_dist)
    print(f"Length of stations in grid: {len(stations)}..\n")
    # plot_stations(eq_lat, eq_lon, stations)

    # sys.exit()
    station_data=[]
    for phase in phases:
        params.phase(phase)
        print(f"Doing phase: {phase}..\n")
        for st in stations:
            params.station(*st)
            try:
                params.geodist('spherical')
                #print(params.asCommandLine(taupserver))
                TimeResult_sph = params.calc(taupserver)
                dist_sph=TimeResult_sph.arrivals[0].distdeg
                time_sph=TimeResult_sph.arrivals[0].time
                params.geodist('geocentric')
                TimeResult_geoc = params.calc(taupserver)
                dist_geoc=TimeResult_geoc.arrivals[0].distdeg
                time_geoc=TimeResult_geoc.arrivals[0].time
                params.geodist('geodetic')
                TimeResult_geod = params.calc(taupserver)
                dist_geod=TimeResult_geod.arrivals[0].distdeg
                time_geod=TimeResult_geod.arrivals[0].time

                station_data.append({
                    'lat': st[0],
                    'lon': st[1],
                    'dist_sph': dist_sph,
                    'time_sph': time_sph,
                    'dist_geod': dist_geod,
                    'time_geod': time_geod,
                    'dist_geoc': dist_geoc,
                    'time_geoc': time_geoc,
                    'phase': phase})
            except:
                print(f'No arrivals for station: {st} {phase} ')

    print(f"lenth of station data: {len(station_data)}")
    return station_data
#


#
def fig6(taupserver):
    spacing=5
    eq_lat=45
    eq_lon=0
    station_data=cal_dist_time(taupserver,phase=['SS'],spacing=spacing, max_dist=175,eq_lon=eq_lon,eq_lat=eq_lat)

    station_data_epy,arr=ellipticipy_py_corr(phase=['SS'],spacing=spacing, max_dist=175,eq_lon=eq_lon,eq_lat=eq_lat)

    ## Fig 6
    # distance differences
    plot_dist_diff(station_data,'dist_geoc','dist_geod',eq_lon=eq_lon,eq_lat=eq_lat,cmap='PRGn',label="$ \Delta $ Distance ($^\\circ$)",tick_space=.1,figname='geod_geoc_dist_d.png')
    plot_dist_diff(station_data,'dist_geoc','dist_sph',eq_lon=eq_lon,eq_lat=eq_lat,cmap='PRGn',label="$ \Delta $ Distance ($^\\circ$)",tick_space=.1,figname='geoc_dist_d.png')
    # ##### time differences
    plot_dist_diff(station_data,'time_geoc','time_geod',eq_lon=eq_lon,eq_lat=eq_lat,cmap='RdGy',label="$ \Delta $ Time (s)",tick_space=.5,figname='geod_geoc_time_d.png')
    plot_dist_diff(station_data,'time_geoc','time_sph',eq_lon=eq_lon,eq_lat=eq_lat,cmap='RdGy',label="$ \Delta $ Time (s)",tick_space=1,figname='geoc_time_d.png')

    # Fig S2.4, supplementary
    ###### time correction using ellipticipy
    plot_dist_diff(station_data_epy,'time_correction','None',eq_lon=eq_lon,eq_lat=eq_lat,cmap='RdGy',label="Ellipticity Correction (s)",tick_space=1,figname='ellip_py_sphr_dist_corr.png')

def main():
    taup_path="~/Research/sct_wat/TauP/build/install/TauP/bin/taup"
    taup_path="~/Code/seis/TauP/build/install/TauP/bin/taup"

    with taup.TauPServer(taup_path=taup_path) as taupserver:
        fig6(taupserver)

if __name__ == '__main__':
    main()
