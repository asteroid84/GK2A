#!/home/hyungmin/anaconda3/bin/python
# file : conv_gk2a_l1b_nc.py

import numpy as np
from netCDF4 import Dataset


## NC File
#ncfile = '/Data04/GK2A_L1B/FD/201909/01/00/gk2a_ami_le1b_wv063_fd020ge_201909010000.nc'

## Constants
NaN = float('nan')
planck_c = 6.62607004E-34  # Planck's constant [J*s]
boltz_c = 1.38064852E-23   # Boltzmann constant [J/K]
light_speed = 2.99792458E+8 # Speed of light in vacuum [m/s]
c1 = 2 * planck_c * (light_speed) ** 2
c2 = (planck_c * light_speed) / boltz_c

## Calibration Table v3.0
lamda_vis = [0.47020000, 0.50860000, 0.63940000, \
             0.86300000, 1.37400000, 1.60920000]
gain_vis  = [0.363545805215835, 0.343625485897064, 0.154856294393539, 
             0.0457241721451282, 0.0346878096461296, 0.0498007982969284]
offset_vis = [-7.27090454101562, -6.87249755859375, -6.19424438476562, \
              -3.65792846679687, -1.38751220703125, -0.996017456054687]
r_to_a_vis = [0.0015582450, 0.0016595767, 0.0019244840, \
              0.0032723873, 0.0087081313, 0.0129512876]

# [um]
lamda_ir = [3.82749133182218, 6.18196269968843, 6.93685519362696, \
            7.32466585491178, 8.58406387549677, 9.61575244772356, \
            10.35032342314430, 11.21436982235590, 12.33640374446540, \
            13.26980420792720]

# wavenumber [cm^-1]
wavenum_ir = [2612.677373521110, 1617.609242531340, 1441.575428760170, \
              1365.249992024440, 1164.949392856340, 1039.960216776110, \
              966.153383926055, 891.713057301260, 810.609007871230, \
              753.590621482278]

gain_ir = [-0.00108296517282724000, -0.01089146733283990000, -0.00818779878318309000, \
           -0.00969827175140380000, -0.01448065508157010000, -0.01784354634582990000, \
           -0.01981969550251960000, -0.02167448587715620000, -0.02337997220456600000, \
           -0.02430375665426250000]

offset_ir = [17.69998741149900000000, 44.17770385742180000000, 66.74807739257810000000, \
             79.06085205078120000000, 118.05090332031200000000, 145.46487426757800000000, \
             161.58013916015600000000, 176.71343994140600000000, 190.64962768554600000000, \
             198.22436523437500000000]

coef_c0 = [-0.447843939824124, -1.762794940111470, -0.334311414359106, \
           -0.061312485969660, -0.141418528203155, -0.114017728158198, \
           -0.142866448475177, -0.249111718496148, -0.458113885722738, \
           -0.093852156852766]

coef_c1 = [1.000655680903890, 1.004149105622780, 1.000973598744680, \
           1.000190087229410, 1.000522329068850, 1.000473805854020, \
           1.000640695720490, 1.001211668737560, 1.002455209755350, \
           1.000539821129660]

coef_c2 = [-6.338240899124480E-08, -9.833109143193850E-07, -4.946030702523040E-07, \
           -1.058636567504990E-07, -3.628727607610900E-07, -3.749315099284030E-07, \
           -5.504432949604980E-07, -1.131679640116650E-06, -2.530643147204760E-06, \
           -5.949137153128490E-07]

min_bound = [20, 20, 40, 80, 40, \
             21, 0, 0, 0, 0, \
             0, 0, 0, 0, 0, \
             0]

max_bound = [2047, 2047, 4095, 8191, 4095, \
#            11bit 11bit 12bit 13bit 12bit
             2047, 16344, 4095, 8152, 8152, \
#            11bit 14bit  12bit 13bit 13bit     
             8152, 8152, 8152, 8153, 8154, 8156]   
#            13bit 13bit 13bit 13bit 13bit 13bit



##=============================================================================
def read_gk2a_l1b_nc(ncfile): 
  
    ch_names = ['vi004', 'vi005', 'vi006', 'vi008', 'nr013', \
                'nr016', 'sw038', 'wv063', 'wv069', 'wv073', \
                'ir087', 'ir096', 'ir105', 'ir112', 'ir123', 'ir133']

    file_name = ncfile.split('/')
    file_name = file_name[-1]
    ch = file_name.split('_')
    ch = ch[3]
    ich = ch_names.index(ch)

    print('GK2A AMI file open: ', ncfile)
    print('           Channel: ', ch)
    print()

    nc_obj = Dataset(ncfile, mode='r')
    var_name = 'image_pixel_values'  
#  dim_x, dim_y = nc_obj.variables[var_name].shape

    count = nc_obj.variables[var_name][:,:]
    bad = np.where((count < min_bound[ich]) | (count > max_bound[ich]))

    count = np.float_(count[:,:])
    count[bad] = float('nan')

    print('Count Boundary at this channel: ', min_bound[ich], ' to ', max_bound[ich])
    print('                   Count Range: ', np.nanmin(count), ' to ', np.nanmax(count))
    print()

    return count, ich
##=============================================================================


##=============================================================================
def convert_gk2a_count2rad(count, ich):
  
    if(ich < 6): # for VIS and N-IR channels (ch 01 ~ 06)
        i = ich
        rad = gain_vis[i] * count + offset_vis[i]

    else:        # for IR channels (ch 07 ~ 16)
        i = ich - 6
        rad = gain_ir[i] * count + offset_ir[i]

    print('                Radiance Range: ', np.nanmin(rad), ' to ', np.nanmax(rad))
    print()

    return rad
##=============================================================================


##=============================================================================
def convert_gk2a_rad2alb_or_tbb(rad, ich):
   
    if(ich < 6): # for VIS and N-IR channels (ch 01 ~ 06)
        alb = rad * r_to_a_vis[ich] * 100
        print('             Albedo Range: ', np.nanmin(alb), ' to ', np.nanmax(alb))

    else:        # for IR channels (ch 07 ~ 16)
        i = ich - 6
        real_temp_eff = (c2 * (wavenum_ir[i] * 100.) / (np.log( (c1 * (wavenum_ir[i] * 100.)**3) / (rad * 1.E-5) + 1.)) )
 
        tbb = (coef_c0[i] + coef_c1[i] * real_temp_eff + coef_c2[i] * real_temp_eff**2)

        print('                    Tbb Range: ', np.nanmin(tbb), ' to ', np.nanmax(tbb))  

    print()
    return tbb
##=============================================================================


##=============================================================================
def convert_gk2a_nc2alb_or_tbb(ncfile):

    count, ich = read_gk2a_l1b_nc(ncfile)
    rad = convert_gk2a_count2rad(count, ich)
    alb_tbb = convert_gk2a_rad2alb_or_tbb(rad, ich)
 
    return alb_tbb
##=============================================================================
  
#count, ich = read_gk2a_l1b_nc(ncfile)
#rad = convert_gk2a_count2rad(count, ich)
#tbb = convert_gk2a_rad2alb_or_tbb(rad, ich)

#alb_tbb = convert_gk2a_nc2alb_or_tbb(ncfile)
