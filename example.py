#!/home/hyungmin/anaconda3/bin/python
# file : example.py

import gk2a  
import numpy as np

ncfile = '/Data04/GK2A_L1B/FD/201909/01/00/gk2a_ami_le1b_wv063_fd020ge_201909010000.nc'
lon_file_2km = '/home/hyungmin/mywork/1_data/3_gk2a/GK2A_AUX/LatLon/Lon_2km.bin'
lat_file_2km = '/home/hyungmin/mywork/1_data/3_gk2a/GK2A_AUX/LatLon/Lat_2km.bin'
lon_file_1km = '/home/hyungmin/mywork/1_data/3_gk2a/GK2A_AUX/LatLon/Lon_1km.bin'
lat_file_1km = '/home/hyungmin/mywork/1_data/3_gk2a/GK2A_AUX/LatLon/Lat_1km.bin'
lon_file_500 = '/home/hyungmin/mywork/1_data/3_gk2a/GK2A_AUX/LatLon/Lon_500m.bin'
lat_file_500 = '/home/hyungmin/mywork/1_data/3_gk2a/GK2A_AUX/LatLon/Lat_500m.bin'


#tbb63 = gk2a.level1b.convert_nc2alb_or_tbb(ncfile)
#rad63 = gk2a.level1b.convert_nc2rad(ncfile)
#lon2km, lat2km = gk2a.read_lonlat.fulldisk(lon_file_2km, lat_file_2km)
#lon1km, lat1km = gk2a.read_lonlat.fulldisk(lon_file_1km, lat_file_1km)
#lon500, lat500 = gk2a.read_lonlat.fulldisk(lon_file_500, lat_file_500)

lon, lat = gk2a.calc_lonlat.colline_to_lonlat('2km', 1000., 1000.)
col, lin = gk2a.calc_lonlat.lonlat_to_colline('2km', lon, lat)

print(lon, lat)
print(col, lin)
