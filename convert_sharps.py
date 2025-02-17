#Script to take downloaded sharp data and output as a magnetic field at a sensible resolution
#Will require some degree of smoothing to make it work reasonably.
#Hopefully can integrate all together in a mf_smart equivalent.

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import netcdf_file
import os
import time

from obtain_sharps import obtain_sharp
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RegularGridInterpolator

def sharp_info(sharp_id, root_fname):
    #For a series of SHARP data, returns the first and last good file, and the resolutions
    raw_directory = root_fname + '%05d_raw/' % sharp_id
    #Check length of this list against info

    flist = os.listdir(raw_directory)
    flist.sort()

    missing_files = []
    #Check all are downloaded
    count = 0
    for fi, fname in enumerate(flist):
        if count != int(fname[-8:-3]):
            missing_fname = raw_directory + '%05d_%05d.nc' % (sharp_id, count)
            print('File in sequence not found:', missing_fname)
            missing_files.append(missing_fname)
            count += 1
        count += 1

    if len(missing_files) > 0:
        print('Attempting to download missing files...')
        obtain_sharp(sharp_id, root_fname, files = missing_files)
        print('Hopefully sucessful. Carrying on.')

    for fi, fname in enumerate(flist):
        try:
            data = netcdf_file(raw_directory + fname, 'r', mmap=False)
            bp = np.swapaxes(data.variables['Bp'][:], 0,1)
            data.close()
            if np.sum(~np.isnan(bp)) == bp.shape[0]*bp.shape[1]:
                if fi == 0:
                    start = fi   #Start not near a limb
                else:
                    start = fi + 80     #Start near a limb
                break
        except:
            pass

    for fi, fname in enumerate(flist[::-1]):
        try:
            data = netcdf_file(raw_directory + fname, 'r', mmap=False)
            bp = np.swapaxes(data.variables['Bp'][:], 0,1)
            data.close()
            if np.sum(~np.isnan(bp)) == bp.shape[0]*bp.shape[1]:
                end = len(flist) - 2 - fi
                if fi == 0:
                    end = len(flist) - 2 - fi
                else:
                    end = len(flist) - 2 - fi - 80
                break
        except:
            pass

    return bp.shape[0], bp.shape[1], start, end

def balance_flux(field):
    #While retaining the sign of each cell, enforces proper flux balance. Can do this in the convert stage, alternatively.
    #This shouldn't be necessary if the initial data has been balanced
    #plt.pcolormesh(field.T, cmap = 'seismic', vmin = -np.max(np.abs(field)), vmax = np.max(np.abs(field)))
    #plt.show()
    # toplot = field
    # plt.pcolormesh(toplot.T, cmap = 'seismic', vmin = -np.max(np.abs(toplot)), vmax = np.max(np.abs(toplot)))
    # plt.show()

    #Attempt to find a way to balance EACH column and row here.

    fieldplus = field[1:-1,1:-1] * [field[1:-1,1:-1] >= 0]
    fieldminus = field[1:-1,1:-1] * [field[1:-1,1:-1] < 0]

    avg = (np.sum(fieldplus) - np.sum(fieldminus))/2

    fieldplus = fieldplus[0] * avg/np.sum(fieldplus)
    fieldminus = fieldminus[0] * -avg/np.sum(fieldminus)
    field = np.array(fieldplus + fieldminus)

    # toplot = field
    # plt.pcolormesh(toplot.T, cmap = 'seismic', vmin = -np.max(np.abs(toplot)), vmax = np.max(np.abs(toplot)))
    # plt.show()

    return np.array(fieldplus + fieldminus)

#Want to import this as a function so just do it like that
def convert_sharp(grid, sharp_id, root_fname, start = 0, end = 1, max_mags = 10000, plot = False, normalise = False, envelope_factor = -1, padding_factor = 0):
    #Imports the grid data, sharp id and the start and end frames to plot
    #Generally will be a mess at the start and end, so ignore those bits.
    print('Converting raw data from SHARP', sharp_id)
    #for import_number in range(start, end):
    output_dir = root_fname + '%05d_mag/' % sharp_id
    raw_directory = root_fname + '%05d_raw/' % sharp_id

    if os.path.exists(output_dir):
        for fname in os.listdir(output_dir):
            os.system('rm ' + output_dir+fname)
    else:
        os.mkdir(output_dir)

    #If there are too many input files, thin them out a bit at this point
    #File admin done. Read in files and do some converting...

    output_count = 0
    thin_fact = max(1, int((end - start)/max_mags) + 1)
    raw_mag_times = np.load('./parameters/raw_times%05d.npy' % sharp_id)

    mag_times = []
    for fi in range(start, end, thin_fact):

        mag_times.append(raw_mag_times[fi])

        fname = raw_directory + '%05d_%05d.nc' % (sharp_id, fi)

        data = netcdf_file(fname, 'r', mmap=False)
        bp = np.swapaxes(data.variables['Bp'][:], 0,1)
        bt = np.swapaxes(data.variables['Bt'][:], 0,1)
        br = np.swapaxes(data.variables['Br'][:], 0,1)
        try:
            data = netcdf_file(fname, 'r', mmap=False)
            bp = np.swapaxes(data.variables['Bp'][:], 0,1)
            bt = np.swapaxes(data.variables['Bt'][:], 0,1)
            br = np.swapaxes(data.variables['Br'][:], 0,1)
            data.close()

        except:
            print('Data error with fname', fname, '. Attempting to re-download...')
            obtain_sharp(sharp_id, root_fname, files = [fname])

            data = netcdf_file(fname, 'r', mmap=False)
            bp = np.swapaxes(data.variables['Bp'][:], 0,1)
            bt = np.swapaxes(data.variables['Bt'][:], 0,1)
            br = np.swapaxes(data.variables['Br'][:], 0,1)
            data.close()

            continue

        #Need to flip the coordinates of bt and bp (just negative) and do some blurring.
        #Amount of blurring still to be decided, but shouldn't matter too much
        grid_downscale = bp.shape[0]/grid.nx

        if grid_downscale < 2.0:
            raise Exception('Interpolating to a very fine grid. Consider not doing that...')

        #NOT SURE WHETHER TO FLIP SIGNS HERE. I THINK IT'S JUST y/theta
        bx_smooth = gaussian_filter(bp, sigma = grid_downscale*0.5)
        by_smooth = gaussian_filter(-bt, sigma = grid_downscale*0.5)
        bz_smooth = gaussian_filter(br, sigma = grid_downscale*0.5)

        #Interpolate onto the STAGGERED grid (just linearly)
        #If padding, need to be more smart here

        pad_distance =  max(grid.x1, grid.y1) - max(grid.x1, grid.y1)/(1.0 + padding_factor)
        #envelope_factor = 1.0/padding_factor

        xs_import = np.linspace(grid.x0+pad_distance, grid.x1-pad_distance, bx_smooth.shape[0])
        ys_import = np.linspace(grid.y0+pad_distance, grid.y1-pad_distance, by_smooth.shape[1])

        X, Y = np.meshgrid(grid.xs, grid.yc, indexing = 'ij')
        bx_fn = RegularGridInterpolator((xs_import, ys_import), bx_smooth, bounds_error = False, method = 'linear', fill_value = 0.0)
        bx_out = bx_fn((X,Y))   #Difference now interpolated to the new grid

        X, Y = np.meshgrid(grid.xc, grid.ys, indexing = 'ij')
        by_fn = RegularGridInterpolator((xs_import, ys_import), by_smooth, bounds_error = False, method = 'linear', fill_value = 0.0)
        by_out = by_fn((X,Y))   #Difference now interpolated to the new grid

        X, Y = np.meshgrid(grid.xc, grid.yc, indexing = 'ij')
        bz_fn = RegularGridInterpolator((xs_import, ys_import), bz_smooth, bounds_error = False, method = 'linear', fill_value = 0.0)
        bz_out = bz_fn((X,Y))   #Difference now interpolated to the new grid

        bx_out = gaussian_filter(bx_out, sigma = 0.5)
        by_out = gaussian_filter(by_out, sigma = 0.5)
        bz_out = gaussian_filter(bz_out, sigma = 0.5)


        # toplot = bz_out.T
        # plt.imshow(toplot,cmap ='seismic', vmax = np.max(np.abs(toplot)), vmin = -np.max(np.abs(toplot)))
        # plt.show()

        if normalise:
            bz_out[1:-1,1:-1] = balance_flux(bz_out)
            if fi == start:
                norm_factor = np.max(np.abs(bz_out))
            bx_out = bx_out/norm_factor
            by_out = by_out/norm_factor
            bz_out = bz_out/norm_factor

        if plot:
            fig, axs = plt.subplots(3,3, figsize = (10,10))

            toplot = bx_smooth.T
            axs[1,0].imshow(toplot,cmap ='seismic', vmax = np.max(np.abs(toplot)), vmin = -np.max(np.abs(toplot)))
            toplot = by_smooth.T
            axs[1,1].imshow(toplot,cmap ='seismic', vmax = np.max(np.abs(toplot)), vmin = -np.max(np.abs(toplot)))
            toplot = bz_smooth.T
            axs[1,2].imshow(toplot,cmap ='seismic', vmax = np.max(np.abs(toplot)), vmin = -np.max(np.abs(toplot)))

            toplot = bp.T
            axs[0,0].imshow(toplot,cmap ='seismic', vmax = np.max(np.abs(toplot)), vmin = -np.max(np.abs(toplot)))
            toplot = bt.T

            axs[0,1].imshow(toplot,cmap ='seismic', vmax = np.max(np.abs(toplot)), vmin = -np.max(np.abs(toplot)))
            toplot = br.T
            axs[0,2].imshow(toplot,cmap ='seismic', vmax = np.max(np.abs(toplot)), vmin = -np.max(np.abs(toplot)))

            toplot = bx_out.T
            axs[2,0].imshow(toplot,cmap ='seismic', vmax = np.max(np.abs(toplot)), vmin = -np.max(np.abs(toplot)))
            toplot = by_out.T
            axs[2,1].imshow(toplot,cmap ='seismic', vmax = np.max(np.abs(toplot)), vmin = -np.max(np.abs(toplot)))
            toplot = bz_out.T
            axs[2,2].imshow(toplot,cmap ='seismic', vmax = np.max(np.abs(toplot)), vmin = -np.max(np.abs(toplot)))

            plt.savefig('./import_plots/%05d' % output_count)
            plt.close()

        #Save out to new filename

        fid = netcdf_file(output_dir + '%04d.nc' % output_count , 'w')
        fid.createDimension('xs', grid.nx+1)
        fid.createDimension('ys', grid.ny+1)
        fid.createDimension('xc', grid.nx+2)
        fid.createDimension('yc', grid.ny+2)

        vid = fid.createVariable('xs', 'd', ('xs',))
        vid[:] = grid.xs
        vid = fid.createVariable('ys', 'd', ('ys',))
        vid[:] = grid.ys

        vid = fid.createVariable('xc', 'd', ('xc',))
        vid[:] = grid.xc
        vid = fid.createVariable('yc', 'd', ('yc',))
        vid[:] = grid.yc

        #Transposes are necessary as it's easier to flip here than in Fortran

        vid = fid.createVariable('bx', 'd', ('yc','xs'))
        vid[:] = np.swapaxes(bx_out, 0, 1)
        vid = fid.createVariable('by', 'd', ('ys','xc'))
        vid[:] = np.swapaxes(by_out, 0, 1)
        vid = fid.createVariable('bz', 'd', ('yc','xc'))
        vid[:] = np.swapaxes(bz_out, 0, 1)
        fid.close()

        output_count += 1

    np.save('./parameters/mag_times%05d.npy' % sharp_id, np.array(mag_times))







