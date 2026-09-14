#downloads and saves Fiji eq data as sac files. Also saves the inventory!
from obspy import read, Stream, UTCDateTime,read_events
from obspy.core.event import Origin, Catalog
from obspy.core.inventory.inventory import read_inventory
import numpy as np
from obspy.clients.fdsn import Client
from obspy.clients.fdsn import RoutingClient
import datetime
import os
####

inv_fiji = "inventory_big_fiji.xml"

def getData():
    iris = RoutingClient("earthscope-federator")
    eida = RoutingClient("eida-routing")

    starttime= UTCDateTime('2018-08-19T00:19:40') #Fiji

    eq_lat,eq_long,eq_depth=-18.1125,-178.1530,600 #Fiji
    endtime= starttime+7000
    network_eu="NS,HE,KO,FR,GR,IV"
    network_list="II,IM,IU,CU,IC,GT,AK,CN,US,AU,GB,GE"#

    try:
        inventory_big=read_inventory(inv_fiji)
    except err:
        print(err)
        print("reloading inventory...")
        inventory_big = iris.get_stations(network=network_list,starttime=starttime,endtime=endtime)
        for net in ["NS", "HE", "KO", "FR", "GR", "IV"]:
            try:
                inventory_big += eida.get_stations(network=net,starttime=starttime,endtime=endtime)
                print(net, "Inventory OK")
            except Exception:
                print(net, "NO DATA")
        inventory_big.write(inv_fiji, format="STATIONXML")
        print(f"inventory saved to {inv_fiji}")

    # inventory_big.plot(label=False,color_per_network=True,resolution='i',continent_fill_color='honeydew',alpha=.5)
    # download data bit
    stream_all = iris.get_waveforms(network=network_list, station='*', location="00",channel= "B*Z", starttime=starttime,endtime=endtime)
    for net in ["NS", "HE", "KO", "FR", "GR",'IV']:
        try:
            stream_all += eida.get_waveforms(network=net, station='*', location="*", channel="B*Z",starttime= starttime,endtime=endtime)#,attach_response=False)
            print(net, "Waveform OK")
        except Exception:
            print(net, "NO DATA")
    print('len of stream:',len(stream_all))
    stream_all.resample(20.0)
    stream_all.filter('bandpass',freqmin=.01, freqmax=.2)

    # print(len(stream_all))
    dirName = "sac_fiji_18"
    os.makedirs(dirName, exist_ok=True)
    for tr in stream_all:
        if tr.stats.npts>130000: # to remove spurious traces with less npts
            sst='{}.sac'.format(tr.id)
            data = os.path.join(dirName, sst)
            tr.write(data,format='sac')
    ####

if __name__ == "__main__":
    getData()
