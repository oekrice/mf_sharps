#An attempt to see if I can download SHARPS data for use with the helicity matching scheme

import numpy as np
import matplotlib.pyplot as plt
import drms
from astropy.io import fits
from scipy.io import netcdf_file
import os

def obtain_sharp(sharpnum, fname_root, plot_flag = 0, files = []):
    series = 'hmi.sharp_cea_720s'
    segments = ['Bp','Bt','Br']
    kwlist = ['T_REC']

    email = 'oliver.e.rice@durham.ac.uk'

    input_folder = './Data/'

    c = drms.Client(email=email) ##Use your own email address.
    k = c.query('%s[%d]' % (series, sharpnum), seg = segments)   #What does this do?!

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

    print(files)
    for fname in files:
        #Save out as netcdf
        if plot_flag:
            fig, axs = plt.subplots(3)

        fid = netcdf_file(fname, 'w')
        fi = int(int(fname[-8:-3]))
        for i in range(3):   #Each component
            url_end = k[segments[i]][fi]
            url = 'http://jsoc.stanford.edu' + url_end
            data = fits.getdata(url)

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

        if plot_flag:

            ax = axs[i]
            ax.imshow(data)

            fid.close()

            plt.savefig('/extra/tmp/trcn27/sharps/plots/%05d_%05d.png' % (sharpnum, fi))

            plt.close()
