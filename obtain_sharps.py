#An attempt to see if I can download SHARPS data for use with the helicity matching scheme

import numpy as np
import matplotlib.pyplot as plt
import drms
from astropy.io import fits
from scipy.io import netcdf_file
import os
from datetime import datetime, timedelta

def lookup_sharp(region_id):
    print('Finding relevant SHARP data...')
    email = 'oliver.e.rice@durham.ac.uk'

    c = drms.Client(email=email) ##Use your own email address.
    k = c.query('hmi.Mharp_720s[][][? NOAA_ARS ~ "%d" ?]' % region_id , key= ['HARPNUM'], n = 1)

    return int(k['HARPNUM'][0])

def time_interval(string1, string2):
    #With the silly time format, gives the number of minutes between each one

    string1_proper = string1[0:4] + string1[5:7] + string1[8:10] + 'T' + string1[11:19]
    string2_proper = string2[0:4] + string2[5:7] + string2[8:10] + 'T' + string2[11:19]

    date1 = datetime.fromisoformat(string1_proper)
    date2 = datetime.fromisoformat(string2_proper)

    return (date2 - date1).days*(60*24) + (date2 - date1).seconds/60


def obtain_sharp(sharpnum, fname_root, plot_flag = 0, files = []):

    series = 'hmi.sharp_cea_720s'
    segments = ['Bp','Bt','Br']
    kwlist = ['T_REC']

    email = 'oliver.e.rice@durham.ac.uk'

    c = drms.Client(email=email) #Use your own email address.
    times, k = c.query('%s[%d]' % (series, sharpnum), seg = segments, key = kwlist)   #What does this do?!

    #Save out the times at this point, just in case. In minutes after the first one.
    mag_times = []

    for ti in range(len(times)):
         if ti == 0:
             mag_times.append(0.)
         else:
             mag_times.append(time_interval(times['T_REC'][0], times['T_REC'][ti]))

    np.save('./parameters/raw_times%05d.npy' % sharpnum, np.array(mag_times))
    allfiles = k[segments[0]]

    nfiles = len(allfiles)

    if len(files) == 0:
        print(nfiles, 'files found to download. Doing so...')

    fname_root = fname_root + '%05d_raw/' % sharpnum

    if not os.path.exists(fname_root):
        os.mkdir(fname_root)

    if len(files) == 0:  #Download everyting
        files = []
        for fi in range(0,nfiles):
            fname = fname_root + '%05d_%05d.nc' % (sharpnum, fi)
            files.append(fname)

    for fi, fname in enumerate(files):

        #Check if file exists, and skip if so. If data is corrupted will sort out later
        if os.path.exists(fname):
            pass

        #Save out as netcdf

        fid = netcdf_file(fname, 'w')
        fi = int(int(fname[-8:-3]))
        for i in range(3):   #Each component
            url_end = k[segments[i]][fi]
            url = 'http://jsoc.stanford.edu' + url_end
            try:
                data = fits.getdata(url)
            except:
                print('Data not found... will try later')
                continue

            nx = data.shape[0]; ny = data.shape[1]
            if i == 0:

                fid.createDimension('xs', nx)
                fid.createDimension('ys', ny)

                vid = fid.createVariable('xs', 'd', ('xs',))
                vid[:] = np.linspace(0,1,nx)
                vid = fid.createVariable('ys', 'd', ('ys',))
                vid[:] = np.linspace(0,1,ny)

            vid = fid.createVariable(segments[i], 'd', ('xs','ys'))
            vid[:] = data

        fid.close()
        pc = fi/len(files)
        print('Downloaded file', fname, '(', int(pc*100), '%)', end='\r')
