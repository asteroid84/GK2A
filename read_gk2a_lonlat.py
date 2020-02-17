#!/home/hyungmin/anaconda3/bin/python
# file : read_gk2a_lonlat.py

import struct
import numpy as np

#lon_file = '/home/hyungmin/mywork/1_data/3_gk2a/GK2A_AUX/LatLon/Lon_2km.bin'
#lat_file = '/home/hyungmin/mywork/1_data/3_gk2a/GK2A_AUX/LatLon/Lat_2km.bin'

def read_gk2a_fulldisk(lon_file, lat_file):
  
    xs = 5500
    ys = 5500

    with open(lon_file,'rb') as binary:
        buf = np.fromfile(binary, dtype='f', count=xs*ys)
        lon = np.reshape(buf, [xs, ys], order='F')

    with open(lat_file,'rb') as binary:
        buf = np.fromfile(binary, dtype='f', count=xs*ys)
        lat = np.reshape(buf, [xs, ys], order='F')

    bad = np.where((lon < -180) | (lon > 180) | (lat < -90) | (lat > 90))
    lon[bad] = float('nan')
    lat[bad] = float('nan')

    print('  GK2A Longitude range: ',np.nanmin(lon),' to ',np.nanmax(lon))
    print('  GK2A  Latitude range: ',np.nanmin(lat),' to ',np.nanmax(lat))    
    print()

    return lon, lat

#lon, lat = read_gk2a_fulldisk(lon_file, lat_file)
