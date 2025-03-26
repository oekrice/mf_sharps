#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 15:16:37 2022

@author: trcn27

These diagnostics are used to test for CMEs etc, using the SNAPSHOTS rather than the usual diagnostics saved out.

This is just because the helicity code is tricky to use otherwise.
"""

import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
import matplotlib
from scipy.io import netcdf_file
from scipy.ndimage import gaussian_filter1d
from scipy.fft import ifft2,fft2,fftfreq, fftshift, ifftshift,ifft,fft
from scipy.signal import find_peaks
#from fltrace import trace_fieldlines

#matplotlib.rcParams['text.usetex'] = True

#Doing (basic) diagnostics from the MF emergence runs
#Will try to automate and get parameters etc.

plot_set_number = 0  #Do this now to make things neat. Yes. Good.
runs = [0]
nsets =  len(runs)

cmap = [
    '#1f77b4',
    '#ff7f0e',
    '#2ca02c',
    '#d62728',
    '#9467bd',
    '#8c564b',
    '#e377c2',
    '#7f7f7f',
    '#bcbd22',
    '#17becf'
]
cs = plt.cm.plasma(np.linspace(0.1,0.9,nsets))
machine_flag = 0

def generate_diags():

    for run in runs:

        def norm2d(vec):
            mag = np.linalg.norm(vec)
            if (mag > 0.0):
                v = vec/mag
            else:
                v = np.array([0, 0])
            return np.array([v[0],v[1],0.0])

        def getFrequencyMatrix(ncells,spacing):
            freqlist1da =np.roll(np.linspace(-ncells[0]/2,ncells[0]/2-1,ncells[0]),round(ncells[0]/2))/(ncells[0]*spacing[0])
            freqlist1db =np.roll(np.linspace(-ncells[1]/2,ncells[1]/2-1,ncells[1]),round(ncells[1]/2))/(ncells[1]*spacing[1])
            return np.array([np.array([np.array([2.0*np.pi*freqlist1da[i],2.0*np.pi*freqlist1db[j]]) for j in range(len(freqlist1db))]) for i  in range(len(freqlist1da))]);

        def getFrequencyMatrixVert(ncells,spacing):
            freqlist1 = np.roll(np.linspace(-ncells[2]/2,ncells[2]/2-1,ncells[2]),round(ncells[2]/2))/(ncells[2]*spacing[2])
            return np.array([2.0*np.pi*freqlist1[i] for i  in range(len(freqlist1))])


        try:
            paras = np.loadtxt('./parameters/variables%03d.txt' % run)   #variables numbered based on run number (up to 1000)
            times = np.loadtxt('./parameters/magtimes%03d.txt' % run)
        except:
            print('Run not found')
            continue

        #Import data nicely, grid only needs doing once...

        data_root = '/extra/tmp/mf3d/%03d/' % run

        fig_width = 15#1.0*(513.11743/72)

        nx = int(paras[1])
        ny = int(paras[2])
        nz = int(paras[3])

        nsnaps = int(paras[5])
        ndiags = int(paras[6])

        hflag = int(paras[22])

        tmax = paras[4]

        xs = np.linspace(paras[12],paras[13], nx+1)
        ys = np.linspace(paras[14],paras[15], ny+1)
        zs = np.linspace(paras[16],paras[17], nz+1)

        xc = 0.5*(xs[1:] + xs[:-1])
        yc = 0.5*(ys[1:] + ys[:-1])
        zc = 0.5*(zs[1:] + zs[:-1])

        dx = xs[1] - xs[0]
        dy = ys[1] - ys[0]
        dz = zs[1] - zs[0]

        xc = 0.5*(xs[1:] + xs[:-1])
        yc = 0.5*(ys[1:] + ys[:-1])
        zc = 0.5*(zs[1:] + zs[:-1])

        dx = xs[1] - xs[0]
        dy = ys[1] - ys[0]
        dz = zs[1] - zs[0]

        class Grid():
            def __init__(self):
                self.x0 = xs[0]; self.x1 = xs[-1]
                self.y0 = ys[0]; self.y1 = ys[-1]
                self.z0 = zs[0]; self.z1 = zs[-1]
                self.nx = nx ; self.ny = ny; self.nz = nz

        nsnaps = 1500
        diagnostic_titles = ['Magnetic Energy', 'Magnetic Helicity', 'Avg. Current', 'Max. Current', 'Open Flux', 'Midway Flux', 'Rope Height']
        ndiags = len(diagnostic_titles)

        alldiags = [[] for _ in range(ndiags)]
        alltimes = []

        print('Calculating diagnostics:', diagnostic_titles)
        end_snap = 0
        for snap in range(0,1500):  #Attempt to load in all the magnetic field data
            try:
                data = netcdf_file('%s%04d.nc' % (data_root, snap), 'r', mmap=False)
            except:
                #print('Failed to find data', '%s%04d.nc' % (data_root, snap))
                continue

            bx = np.swapaxes(data.variables['bx'][:],0,2)
            by = np.swapaxes(data.variables['by'][:],0,2)
            bz = np.swapaxes(data.variables['bz'][:],0,2)

            jx = np.swapaxes(data.variables['jx'][:],0,2)
            jy = np.swapaxes(data.variables['jy'][:],0,2)
            jz = np.swapaxes(data.variables['jz'][:],0,2)

            data.close()

            #Calculate diagnostics from these

            bx0 = 0.5*(bx[1:,:,:] + bx[:-1,:,:])
            by0 = 0.5*(by[:,1:,:] + by[:,:-1,:])
            bz0 = 0.5*(bz[:,:,1:] + bz[:,:,:-1])

            jx0 = 0.25*(jx[:,1:,1:] + jx[:,1:,:-1] + jx[:,:-1,1:] + jx[:,:-1,:-1])
            jy0 = 0.25*(jy[1:,:,1:] + jy[1:,:,:-1] + jy[:-1,:,1:] + jy[:-1,:,:-1])
            jz0 = 0.25*(jz[1:,1::,] + jz[1:,:-1:,] + jz[:-1,1::,] + jz[:-1,:-1:,])

            #________________________________
            #TIME
            if snap < len(times):
                alltimes.append(times[snap])
            else:
                alltimes.append(times[-1]*snap/len(times))

            end_snap = snap
            #_________________________________
            #MAGNETIC ENERGY

            if True:
                b0 = bx0**2 + by0**2 + bz0**2
                energy = 0.5*np.sum(b0)*dx*dy*dz
                di = diagnostic_titles.index('Magnetic Energy')
                alldiags[di].append(energy)

            #_________________________________
            #AVG, MAX CURRENT
            if True:
                j0 = jx0**2 + jy0**2 + jz0**2
                sumj = np.sum(j0)*dy*dy*dz
                avgj = np.sum(j0)/(nx*ny*nz)
                di = diagnostic_titles.index('Avg. Current')
                alldiags[di].append(avgj)

                maxj = np.max(j0)
                di = diagnostic_titles.index('Max. Current')
                alldiags[di].append(maxj)

            #_________________________________
            #OPEN FLUX
            if True:
                oflux = np.sum(np.abs(bz[:,:,-1]))
                di = diagnostic_titles.index('Open Flux')
                alldiags[di].append(oflux)

                oflux = np.sum(np.abs(bz[:,:,int(nz//2)]))
                di = diagnostic_titles.index('Midway Flux')
                alldiags[di].append(oflux)

            #________________________________
            #VOLUME HELICITY (IN WINDING GAUGE)
            if True:
                #Using Chris' original code from some flh stuff. Try to use as-is.
                #THIS IS ACTUALLY THE CODE TO GET A FOR THE WINDING
                ncells = [nx, ny, nz]
                spacing = [dx, dy, dz]
                fm = getFrequencyMatrix(ncells,spacing);
                # make the basis
                normFunc = np.vectorize(norm2d)
                kparr = np.array([np.array([norm2d(fm[i][j]) for j in range(len(fm[0]))]) for i  in range(len(fm))]);
                kperp = np.array([np.array([np.array([-kparr[i][j][1],kparr[i][j][0],0.0]) for j in range(len(fm[0]))]) for i  in range(len(fm))])

                # note in the k matrix below the k=0 element is set to one so we can divide by it.
                k = np.array([np.array([1.0 if i==j==0 else np.linalg.norm(fm[i][j]) for i in range(len(fm))]) for j  in range(len(fm[0]))]).T
                Ax = np.zeros(bx0.shape)
                Ay = np.zeros(by0.shape)
                Az = np.zeros(bz0.shape)
                B0fx = np.zeros(bz0.shape[2],dtype=np.complex128)
                B0fy = np.zeros(bz0.shape[2],dtype=np.complex128)
                aftx = np.zeros([bz0.shape[0],bz0.shape[1],bz0.shape[2]],dtype=np.complex128)
                afty  =np.zeros([bz0.shape[0],bz0.shape[1],bz0.shape[2]],dtype=np.complex128)
                aftz = np.zeros([bz0.shape[0],bz0.shape[1],bz0.shape[2]],dtype=np.complex128)
                for i in range(bz0.shape[2]):
                    fbx = fft2(bx0[:,:,i]); fby =fft2(by0[:,:,i]); fbz = fft2(bz0[:,:,i])
                    #set the zero element of the transform for the final vertical part
                    B0fx[i] = fbx[0][0]
                    B0fy[i] = fby[0][0]
                    akperp = -1j*fbz/k
                    ## fix i =j  element
                    akperp[0][0]=0.0
                    akw = 1j*(-(kparr[:,:,1])*fbx + (kparr[:,:,0])*fby)/k
                    ## fix i =j  element
                    akw[0][0]=0.0
                    aftx[:,:,i] = akperp*kperp[:,:,0]
                    afty[:,:,i] = akperp*kperp[:,:,1]
                    aftz[:,:,i] = akperp*kperp[:,:,2]+akw
                # now do the vertical transforms
                kz = getFrequencyMatrixVert(ncells,spacing)
                fB0fx = fft(B0fx)
                fB0fy = fft(B0fy)
                kz[0]=1.0
                ax00 = -1j*fB0fy/kz
                ay00 =  1j*fB0fx/kz
                ax00[0]=0.0
                ay00[0]=0.0
                ax00 = ifft(ax00)
                ay00 = ifft(ay00)
                #finally transform back to real space
                for i in range(bz0.shape[2]):
                    aftx[0][0] = ax00[i]
                    afty[0][0] = ay00[i]
                    ax = ifft2(aftx[:,:,i])
                    ay = ifft2(afty[:,:,i])
                    az = ifft2(aftz[:,:,i])
                    Ax[:,:,i] = np.real(ax)
                    Ay[:,:,i] = np.real(ay)
                    Az[:,:,i] = np.real(az)

                #Sum up the total helicity (A.B), then square root
                helicity = np.abs(Ax*bx0 + Ay*by0 + Az*bz0)
                hsum = np.sqrt(np.sum(helicity*dx*dy*dz))

                di = diagnostic_titles.index('Magnetic Helicity')
                alldiags[di].append(hsum)

            #________________________________
            #ROPE HEIGHT (ONLY WORKS WHEN NICE AND SYMMETRICAL)
            toplot = bx[nx//2, :,:].T
            toplot = bx[:, ny//2,:].T
            if True:
                #If looking at Bx, this will be negative everywhere for an arcade but goes from postive to negative otherwise -- find intercept.
                isrope = False
                di = diagnostic_titles.index('Rope Height')

                bx_slice = 0.5*(bx[nx//2,ny//2,:] + bx[nx//2,ny//2 + 1,:])

                bx_slice[-1] = 0.0

                peaks, _ = find_peaks(-bx_slice, prominence = max(1e-6, 1e-4*np.max(np.abs(bx_slice))))

                if len(peaks) > 0.0:

                    minpoint = sorted(peaks)[-1]

                    if minpoint < 3:
                        isrope = False

                    else:
                        signs = np.array([bx_slice[1:minpoint] > 0], dtype = 'float')[0]*2. - 1.

                        crosses = signs[1:]*signs[:-1]

                        if np.min(crosses) < 0.0:
                            maxcross = np.max(np.where(crosses < 0.0)[0])

                            #Inteprolate between
                            x0 = zc[maxcross + 1]
                            x1 = zc[maxcross + 2]
                            y0 = bx_slice[maxcross + 1]
                            y1 = bx_slice[maxcross + 2]

                            xpt = (x1*abs(y0) + x0*abs(y1))/(abs(y1) + abs(y0))

                            ropeheight = xpt
                            isrope = True
                            plt.scatter(maxcross+1, bx_slice[maxcross+1], c ='red')
                            plt.scatter(peaks, bx_slice[peaks])
                            plt.plot(bx_slice)
                            plt.close()

                        else:
                            isrope = False

                if not isrope:
                    alldiags[di].append(0.0)
                else:
                    alldiags[di].append(ropeheight)
            print('Snap', snap, end='\r')

            np.save('./hdata/diags%d.npy' % run, alldiags)
            np.save('./hdata/times%d.npy' % run, alltimes)
        print('Snap', snap)
        print('Hdata in', alldiags)

    return

def analyse_diags():

    diagnostic_titles = ['Magnetic Energy', 'Magnetic Helicity', 'Avg. Current', 'Max. Current', 'Open Flux', 'Midway Flux', 'Rope Height']

    print('All done, plotting:')

    titles = runs
    cs = ['Red', 'Green', 'Blue', 'Purple']
    for ri, run in enumerate(runs):
        alldiags = np.load('./hdata/diags%d.npy' % run)
        alltimes = np.load('./hdata/times%d.npy' % run)

        print('Hdata out', alldiags)
        ndiags = len(alldiags)

        if ri == 0:

            ncols = int(np.sqrt(ndiags)) + 1
            nrows = int(ndiags/ncols + 0.95)

            fig, axs = plt.subplots(nrows, ncols, figsize = (10, 7))

        if False:
            time_convert = 1.0/(0.05*60)
        else:
            time_convert = 1.0

        for di in range(0,ndiags):

            ax = axs[di//ncols,di%ncols]
            toplot = alldiags[di][:]
            #toplot = gaussian_filter1d(toplot, 2)

            if ri > 0:
                maxy = max(ax.get_ylim()[1], np.percentile(toplot, 99)*1.1)
            else:
                maxy = np.percentile(toplot, 99)*1.1
            ax.set_ylim(0.0, maxy)
            ax.plot(alltimes[:]*time_convert, toplot, c = cmap[ri])

            ax.set_title(diagnostic_titles[di])

        ax = axs[-1,-1]
        ax.axhline(0.0, label = titles[ri], c = cmap[ri])
        plt.legend()

    plt.tight_layout()
    plt.savefig('diagnostics.png')
    plt.show()

generate_diags()
analyse_diags()






