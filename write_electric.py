#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 15:16:37 2022

@author: trcn27

Generates .nc files of the electric field, to be read-in to the Fortran code.
Needs to deal with different resolutions as the inputs are all the same (192 I think)
Read in the 'magnetograms' from the .nc files

THIS VERSION DIFFERS:
No interpolation required as it uses sharp data --

"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib
from scipy.io import netcdf_file
from scipy.interpolate import RegularGridInterpolator
from scipy.fft import fft, ifft2, fft2, ifft
import os
from scipy.ndimage import gaussian_filter
from scipy.linalg import eigh_tridiagonal

def read_boundary(fname):
    #Reads in saved magnetic field and just outputs bx, by, bz
    data = netcdf_file(fname, 'r', mmap=False)
    #Check if this is just the boundary data or everything
    if len(data.variables['bx'][:].shape) == 3:
        bx = np.swapaxes(data.variables['bx'][:], 0, 2)[1:-1,:,0]
        by = np.swapaxes(data.variables['by'][:], 0, 2)[:,1:-1,0]
        bz = np.swapaxes(data.variables['bz'][:], 0, 2)[:,:,0]
    else:
        bx = np.swapaxes(data.variables['bx'][:], 0, 1)[1:-1,:]
        by = np.swapaxes(data.variables['by'][:], 0, 1)[:,1:-1]
        bz = np.swapaxes(data.variables['bz'][:], 0, 1)[:,:]
    return bx, by, bz


def compute_inplane_helicity(grid, bx, by, bz, plot = False):
    #Need to average these to grid centres to get the FFT to work
    #Probably should just do on the interior - so resolution at the raw values

    bx0 = 0.5*(bx[1:,1:-1] + bx[:-1,1:-1])
    by0 = 0.5*(by[1:-1,1:] + by[1:-1,:-1])
    bz0 = bz[1:-1,1:-1]

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

    #def curla(ax, ay, az):
        #Outputs the curl of a (which is averaged to grid points so a bit messy)

    #Find in -plane vector potential in the winding gauge

    fm = getFrequencyMatrix([bz0.shape[0], bz0.shape[1]],[grid.dx, grid.dy]);
    # make the basis

    kparr = np.array([np.array([norm2d(fm[i][j]) for j in range(len(fm[0]))]) for i  in range(len(fm))]);
    kperp = np.array([np.array([np.array([-kparr[i][j][1],kparr[i][j][0],0.0]) for j in range(len(fm[0]))]) for i  in range(len(fm))])
    # note in the k matrix below the k=0 element is set to one so we can divide by it.
    k = np.array([np.array([1.0 if i==j==0 else np.linalg.norm(fm[i][j]) for i in range(len(fm))]) for j  in range(len(fm[0]))]).T

    nx = bz0.shape[0]; ny = bz0.shape[1]
    aftx = np.zeros([bz0.shape[0],bz0.shape[1]],dtype=np.complex128)
    afty = np.zeros([bz0.shape[0],bz0.shape[1]],dtype=np.complex128)
    aftz = np.zeros([bz0.shape[0],bz0.shape[1]],dtype=np.complex128)

    fbx = fft2(bx0[:,:]); fby = fft2(by0[:,:]); fbz = fft2(bz0[:,:])

    akperp = -1j*fbz/k
    ## fix i =j  element
    akw = 1j*(-(kparr[:,:,1])*fbx + (kparr[:,:,0])*fby)/k
    ## fix i =j  element
    aftx[:,:] = akperp*kperp[:,:,0]
    afty[:,:] = akperp*kperp[:,:,1]
    aftz[:,:] = akperp*kperp[:,:,2]+akw

    ax0 = ifft2(aftx[:,:])
    ay0 = ifft2(afty[:,:])
    az0 = ifft2(aftz[:,:])
    ax0 = np.real(ax0)
    ay0 = np.real(ay0)
    az0 = np.real(az0)

    ax = np.zeros((nx, ny+1))
    ay = np.zeros((nx+1, ny))

    ax[:,1:-1] = 0.5*(ax0[:,1:] + ax0[:,:-1])
    ay[1:-1,:] = 0.5*(ay0[1:,:] + ay0[:-1,:])

    bz_test = (ay[1:,:] - ay[:-1,:])/grid.dx - (ax[:,1:] - ax[:,:-1])/grid.dy
    #This vector potential should be reasonably OK... Need code to test though
    if plot:
        fig, axs = plt.subplots(3)
        ax = axs[0]
        im = ax.pcolormesh(bz)
        plt.colorbar(im, ax = ax)

        ax = axs[1]
        im = ax.pcolormesh(bz_test)
        plt.colorbar(im, ax = ax)

        ax = axs[2]
        im = ax.pcolormesh(bz_test[1:-1,1:-1] - bz[2:-2,2:-2])
        plt.colorbar(im, ax = ax)
        plt.tight_layout()
        plt.show()
    #Should be proportional to the magnetic field strength, so this helicity requires a square root. I'm pretty sure the scaling is OK here...
    hfield = np.sqrt(np.abs(ax0*bx0 + ay0*by0 + az0*bz0))

    return hfield


class Grid():
    """In the interest of doing it properly, put grid parameters in here"""
    def __init__(self, run):
        
        paras = np.loadtxt('parameters/variables%03d.txt' % run)

        #Define the grid onto which the electric field should be outputted
        self.nx = int(paras[1])
        self.ny = int(paras[2])

        self.x0 = paras[12]; self.x1 = paras[13]
        self.y0 = paras[14]; self.y1 = paras[15]
                
        self.xs = np.linspace(self.x0, self.x1, self.nx+1)
        self.ys = np.linspace(self.y0, self.y1, self.ny+1)
        
        self.dx = self.xs[1] - self.xs[0]
        self.dy = self.ys[1] - self.ys[0]
        
        self.xc = np.linspace(self.x0 - self.dx/2, self.x1 + self.dx/2, self.nx+2)
        self.yc = np.linspace(self.y0 - self.dy/2, self.y1 + self.dy/2, self.ny+2)

        self.dx = self.xs[1] - self.xs[0]
        self.dy = self.ys[1] - self.ys[0]

    def lap_xribs(self, xribs):
        """Calculates laplacian of a quantity saved on the x ribs, using staggered as appropriate"""
        lap_x = np.zeros((self.nx+2, self.ny+1))
        #x direction
        lap_x[1:-1,1:-1] += (xribs[2:,1:-1] -  2*xribs[1:-1,1:-1] + xribs[:-2,1:-1])/self.dx**2
        lap_x[1:-1,1:-1] += (xribs[1:-1,2:] -  2*xribs[1:-1,1:-1] + xribs[1:-1,:-2])/self.dy**2
        return lap_x

    def lap_yribs(self, yribs):
        """Calculates laplacian of a quantity saved on the y ribs, using staggered as appropriate"""
        lap_y = np.zeros((self.nx+1, self.ny+2))
        #x direction
        lap_y[1:-1,1:-1] += (yribs[2:,1:-1] -  2*yribs[1:-1,1:-1] + yribs[:-2,1:-1])/self.dx**2
        lap_y[1:-1,1:-1] += (yribs[1:-1,2:] -  2*yribs[1:-1,1:-1] + yribs[1:-1,:-2])/self.dy**2
        return lap_y


    def curl_inplane(self, C):
        """Calculates the in-plane curl of the quantity C(x,y)"""
        """Outputs on the respective ribs, hopefully"""
        
        curl_x = np.zeros((self.nx+2, self.ny+1))
        curl_y = np.zeros((self.nx+1, self.ny+2))

        curl_x = (C[:,1:] - C[:,:-1])/self.dy
        curl_y = -(C[1:,:] - C[:-1,:])/self.dx
        
        return curl_x, curl_y

    def curl_E(self, E_x, E_y):
        """Returns the in-plane curl of the vectors E"""
        curl = np.zeros((self.nx+2, self.ny+2))
        curl[1:-1,1:-1] += (E_x[1:-1,1:] - E_x[1:-1,:-1])/self.dy
        curl[1:-1,1:-1] -= (E_y[1:,1:-1] - E_y[:-1,1:-1])/self.dx

        return curl
        
    def lap_points(self, points):
        """Calculates laplacian of a quantity saved on the y ribs, using staggered as appropriate"""
        lap_p = np.zeros((self.nx+1, self.ny+1))
        #x direction
        lap_p[1:-1,1:-1] += (points[2:,1:-1] -  2*points[1:-1,1:-1] + points[:-2,1:-1])/self.dx**2
        lap_p[1:-1,1:-1] += (points[1:-1,2:] -  2*points[1:-1,1:-1] + points[1:-1,:-2])/self.dy**2
        return lap_p

    def lap_centres(self, points):
        """Calculates laplacian of a quantity saved on the y ribs, using staggered as appropriate"""
        lap_p = np.zeros((self.nx+2, self.ny+2))
        #x direction
        lap_p[1:-1,1:-1] += (points[2:,1:-1] -  2*points[1:-1,1:-1] + points[:-2,1:-1])/self.dx**2
        lap_p[1:-1,1:-1] += (points[1:-1,2:] -  2*points[1:-1,1:-1] + points[1:-1,:-2])/self.dy**2
        return lap_p

    def div_E(self, E_x, E_y):
        """Returns the in-plane curl of the vectors E"""
        div = np.zeros((self.nx+1, self.ny+1))
        div[:,:] += (E_x[1:,:] - E_x[:-1,:])/self.dx
        div[:,:] += (E_y[:,1:] - E_y[:,:-1])/self.dy

        return div
    
    def grad(self, phi):
        """Returns gradients of phi (at grid points)"""
        grad_x = np.zeros((self.nx+2, self.ny+1))
        grad_y = np.zeros((self.nx+1, self.ny+2))
        grad_x[1:-1,:] = (phi[1:,:] - phi[:-1,:])/self.dx
        grad_y[:,1:-1] = (phi[:,1:] - phi[:,:-1])/self.dy
        
        return grad_x, grad_y
        
        
class FT():
    
    """Creates various quantities related to the Fourier transform (mutliplication matrix etc.)"""
    def __init__(self,grid):
        self.grid = grid
        pass
        
        
    def point_transform(self, rhs_points):
        """2D transform in at the grid points"""
        rhs_transform = ifft2(rhs_points)  
        
        j = np.arange(self.grid.nx+1)
        k = np.arange(self.grid.ny+1)

        d2Xdx2 = (2*np.cos(2*np.pi*j/(self.grid.nx+1)) - 2)/self.grid.dx**2
        d2Ydy2 = (2*np.cos(2*np.pi*k/(self.grid.ny+1)) - 2)/self.grid.dy**2
        
        d2 = d2Xdx2.reshape((self.grid.nx+1,1)) + d2Ydy2
        
        d2[0,0] = 1.0   #Avoid problems with the zero value. All this does is add a constant anyway so set to zero in the transform
        
        a_twiddle = rhs_transform/d2
        
        a_twiddle[0,0] = 0.0
        
        phi = (fft2(a_twiddle)).real
        
        return phi

    def centre_transform(self, rhs_centres):
        """2D transform in at the grid centres"""
        rhs_transform = ifft2(rhs_centres)  
        
        j = np.arange(self.grid.nx+2)
        k = np.arange(self.grid.ny+2)

        d2Xdx2 = (2*np.cos(2*np.pi*j/(self.grid.nx+2)) - 2)/self.grid.dx**2
        d2Ydy2 = (2*np.cos(2*np.pi*k/(self.grid.ny+2)) - 2)/self.grid.dy**2
        
        d2 = d2Xdx2.reshape((self.grid.nx+2,1)) + d2Ydy2
        
        d2[0,0] = 1.0   #Avoid problems with the zero value. All this does is add a constant anyway so set to zero in the transform
        
        a_twiddle = rhs_transform/d2
        
        a_twiddle[0,0] = 0.0
        
        G = (fft2(a_twiddle)).real
        
        return G

def balance_flux(field):
    #While retaining the sign of each cell, enforces proper flux balance. Can do this in the convert stage, alternatively.
    #This shouldn't be necessary if the initial data has been balanced
    fieldplus = field[1:-1,1:-1] * [field[1:-1,1:-1] >= 0]
    fieldminus = field[1:-1,1:-1] * [field[1:-1,1:-1] < 0]

    avg = (np.sum(fieldplus) - np.sum(fieldminus))/2

    fieldplus = fieldplus[0] * avg/np.sum(fieldplus)
    fieldminus = fieldminus[0] * -avg/np.sum(fieldminus)

    return np.array(fieldplus + fieldminus)

def find_eigenthings(grid):
    #Finds the eigenvectors and values for the grid, using the outflow field method
    #x first
    d = -2*np.ones(grid.nx)
    e = np.ones(grid.nx-1)
    #BCs
    d[0] = -1
    d[-1] = -1
    w_x, v_x = eigh_tridiagonal(d,e)
    w_x = w_x[::-1]
    v_x = v_x[:,::-1]

    #Then y
    d = -2*np.ones(grid.ny)
    e = np.ones(grid.ny-1)
    #BCs
    d[0] = -1
    d[-1] = -1
    w_y, v_y = eigh_tridiagonal(d,e)
    w_y = w_y[::-1]
    v_y = v_y[:,::-1]

    return w_x, v_x, w_y, v_y

def transform(rhs, basis, n_x, n_y):
    #Assuming the test vectors are orthogonal, finds the component which matches the rhs
    #Now in 2D...
    num = np.sum(rhs*basis)
    return num

class compute_electrics_bounded():
    #Similarly to the original, but doesn;t use a Fourier transform as the domain has closed boundaries
    def __init__(self, run, init_number, mag_root, mag_times, omega = 0., start = 0, end = 500, initialise = True, plot = False):

        grid = Grid(run)  #Establish grid (on new scales)

        omega_init = omega
        data_directory = mag_root

        paras = np.loadtxt('parameters/variables%03d.txt' % run)

        if not os.path.exists('efields'):
            os.mkdir('efields')

        if not os.path.exists('efields/%03d' % init_number):
            os.mkdir('efields/%03d' % init_number)

        if initialise:   #Need to establish reference helicities from the magnetograms.
            hrefs = []; trefs = []
        elif start == 0:
            halls = []; tplots = []; omegas = []
        else:
            halls = np.load('./hdata/halls%03d.npy' % run).tolist()
            omegas = np.load('./hdata/omegas%03d.npy' % run).tolist()
            tplots = np.load('./hdata/tplots%03d.npy' % run).tolist()

        #Calculate eigenthings (only once)
        w_x, v_x, w_y, v_y = find_eigenthings(grid)
        print('Eigenthings found')

        for snap in range(start, end):

            print('Boundary transform', snap)
            bfield_fname = '%s%04d.nc' % (data_directory, snap)
            efield_fname = '%s%03d/%04d.nc' % ('./efields/', init_number, snap)

            try:
                data = netcdf_file(bfield_fname, 'r', mmap=False)
                #print('File', bfield_fname, 'found')

            except:
                print('File', bfield_fname, 'not found')
                raise Exception('Input Magnetogram Not Found')

                continue

            bfield1 = np.swapaxes(data.variables['bz'][:], 0, 1)
            bfield1[1:-1,1:-1] = balance_flux(bfield1)

            bfield_fname = '%s%04d.nc' % (data_directory, snap + 1)

            try:
                data = netcdf_file(bfield_fname, 'r', mmap=False)

            except:
                print('File', bfield_fname, 'not found')
                raise Exception('Input Magnetogram Not Found')
                continue

            bfield2 = np.swapaxes(data.variables['bz'][:], 0, 1)

            bfield2[1:-1,1:-1] = balance_flux(bfield2)

            if True:
                diff = bfield2 - bfield1   #Difference between the magnetic field at import resolution
            else:
                diff = 0.0*bfield1  #Keep lower boundary constant, for stability testing

            #'Vector potential' is initially called 'G'
            G = np.zeros((grid.nx+2, grid.ny+2))
            #comp = np.zeros((grid.nx+2, grid.ny+2))

            for n_x in range(0,len(w_x)):
                for n_y in range(0,len(w_y)):
                    if n_x + n_y > 0:
                        basis = v_x[:,n_x][:,np.newaxis]*v_y[:,n_y][np.newaxis,:]
                        G[1:-1,1:-1] = G[1:-1,1:-1] + (1.0/(w_x[n_x]/grid.dx**2 + w_y[n_y]/grid.dy**2))*basis*transform(-diff[1:-1,1:-1], basis, n_x, n_y)
                        #comp[1:-1,1:-1] = comp[1:-1,1:-1] + (1.0)*basis*transform(-diff[1:-1,1:-1], basis, n_x, n_y)

            G[0,:] = G[1,:]; G[-1,:] = G[-2,:]
            G[:,0] = G[:,1]; G[:,-1] = G[:,-2]

            lap_test = grid.lap_centres(G)

            ex, ey = grid.curl_inplane(G)
            curl_test = grid.curl_E(ex, ey)

            if plot:
                fig, axs = plt.subplots(3)

                im = axs[0].pcolormesh(curl_test[1:-1,1:-1] )
                plt.colorbar(im, ax = axs[0])
                im = axs[1].pcolormesh(-diff[1:-1,1:-1] )
                plt.colorbar(im, ax = axs[1])
                im = axs[2].pcolormesh(curl_test[1:-1,1:-1] + diff[1:-1,1:-1] )
                plt.colorbar(im, ax = axs[2])

                plt.tight_layout()
                plt.show()
            #Define distribution of the divergence of the electric field.
            #Following Cheung and De Rosa, just proportional to the vertical (averaged to points) magnetic field

            # bf = 0.5 * (bfield1 + bfield2)
            #
            # bf = 0.25*(bf[1:,1:] + bf[:-1,1:] + bf[:-1,:-1] + bf[1:,:-1])
            # D = omega*(bf)
            #
            # div_test = grid.div_E(ex, ey)
            # phi = ft.point_transform(-div_test + D)
            # correct_x, correct_y = grid.grad(phi)
            # ex += correct_x
            # ey += correct_y
            #
            # div_test = grid.div_E(ex, ey)

            if plot:

                if np.max(np.abs(curl_test[1:-1,1:-1] + diff[1:-1,1:-1])) > 1e-10:
                    raise Exception('Electric field calculation failed')
                if np.max(np.abs(div_test[1:-1,1:-1] - D[1:-1,1:-1])) > 1e-10:
                    raise Exception('Electric field calculation failed')

                print(bfield_fname)
                print('Curl test', np.max(np.abs(curl_test[1:-1,1:-1] + diff[1:-1,1:-1])))
                print('Div Test', np.max(np.abs(div_test[1:-1,1:-1] - D[1:-1,1:-1])))
                print('____________________________________________')

            curl_test = grid.curl_E(ex, ey)

            #print('dx', grid.dx, 'dy', grid.dy)
            #Swap sign of E, as that seems to be the way forward.
            ex = -ex
            ey = -ey

            if False:
                plt.pcolormesh(ey, cmap = 'seismic', vmin = -np.max(np.abs(ey)), vmax = np.max(np.abs(ey)))
                plt.savefig('plots/ey%d.png' % snap)
                plt.close()

                plt.pcolormesh(ex, cmap = 'seismic', vmin = -np.max(np.abs(ex)), vmax = np.max(np.abs(ex)))
                plt.savefig('plots/ex%d.png' % snap)
                plt.show()

            fid = netcdf_file(efield_fname, 'w')
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
            #Still doesn't seem to like it -- sums are all correct but going in with the wrong direction.
            vid = fid.createVariable('ex', 'd', ('ys','xc'))
            vid[:] = np.swapaxes(ex, 0, 1)
            vid = fid.createVariable('ey', 'd', ('yc','xs'))
            vid[:] = np.swapaxes(ey, 0, 1)

            fid.close()

            if initialise:  #Calculate reference helicities
                bfield_fname = '%s%04d.nc' % (data_directory, snap)
                data = netcdf_file(bfield_fname, 'r', mmap=False)
                bx = np.swapaxes(data.variables['bx'][:], 0, 1)
                by = np.swapaxes(data.variables['by'][:], 0, 1)
                bz = np.swapaxes(data.variables['bz'][:], 0, 1)
                href = compute_inplane_helicity(grid, bx, by, bz)
                hrefs.append(np.sum(href))
                trefs.append(mag_times[snap])

        if initialise:  #Calculate final one
            bfield_fname = '%s%04d.nc' % (data_directory, snap + 1)
            data = netcdf_file(bfield_fname, 'r', mmap=False)
            bx = np.swapaxes(data.variables['bx'][:], 0, 1)
            by = np.swapaxes(data.variables['by'][:], 0, 1)
            bz = np.swapaxes(data.variables['bz'][:], 0, 1)
            href = compute_inplane_helicity(grid, bx, by, bz)
            hrefs.append(np.sum(href))
            trefs.append(mag_times[snap+1])

        print('Electric fields', start, ' to ', end, ' calculated and saved.')
        if initialise:
            np.save('./hdata/h_ref.npy', np.array(hrefs))
            np.save('./hdata/t_ref.npy', np.array(trefs))
        else:
            self.halls = halls
            self.tplots = tplots
            self.omegas = omegas

class compute_electrics():
    
    def __init__(self, run, init_number, mag_root, mag_times, omega = 0., start = 0, end = 500, initialise = True, plot = False):

        grid = Grid(run)  #Establish grid (on new scales)
        

        omega_init = omega
        data_directory = mag_root
        
        paras = np.loadtxt('parameters/variables%03d.txt' % run)
        
        if not os.path.exists('efields'):
            os.mkdir('efields')

        if not os.path.exists('efields/%03d' % init_number):
            os.mkdir('efields/%03d' % init_number)

        if initialise:   #Need to establish reference helicities from the magnetograms.
            hrefs = []; trefs = []
        elif start == 0:
            halls = []; tplots = []; omegas = []
        else:
            halls = np.load('./hdata/halls%03d.npy' % run).tolist()
            omegas = np.load('./hdata/omegas%03d.npy' % run).tolist()
            tplots = np.load('./hdata/tplots%03d.npy' % run).tolist()

        for snap in range(start, end):

            bfield_fname = '%s%04d.nc' % (data_directory, snap)
            efield_fname = '%s%03d/%04d.nc' % ('./efields/', init_number, snap)
        
            try:
                data = netcdf_file(bfield_fname, 'r', mmap=False)
                #print('File', bfield_fname, 'found')
        
            except:
                print('File', bfield_fname, 'not found')
                raise Exception('Input Magnetogram Not Found')

                continue
        
            bfield1 = np.swapaxes(data.variables['bz'][:], 0, 1)
            bfield1 = balance_flux(bfield1)

            bfield_fname = '%s%04d.nc' % (data_directory, snap + 1)
        
            try:
                data = netcdf_file(bfield_fname, 'r', mmap=False)
        
            except:
                print('File', bfield_fname, 'not found')
                raise Exception('Input Magnetogram Not Found')
                continue
        
            bfield2 = np.swapaxes(data.variables['bz'][:], 0, 1)

            bfield2 = balance_flux(bfield2)

            if True:
                diff = bfield2 - bfield1   #Difference between the magnetic field at import resolution
            else:
                diff = 0.0*bfield1  #Keep lower boundary constant, for stability testing

            ft = FT(grid)
            G = ft.centre_transform(-diff)
            ex, ey = grid.curl_inplane(G)
            curl_test = grid.curl_E(ex, ey)
        
            if plot:
                fig, axs = plt.subplots(3)

                im = axs[0].pcolormesh(curl_test[1:-1,1:-1] )
                plt.colorbar(im, ax = axs[0])
                im = axs[1].pcolormesh(-diff[1:-1,1:-1] )
                plt.colorbar(im, ax = axs[1])
                im = axs[2].pcolormesh(curl_test[1:-1,1:-1] + diff[1:-1,1:-1] )
                plt.colorbar(im, ax = axs[2])

                plt.tight_layout()
                plt.show()
            #Define distribution of the divergence of the electric field. 
            #Following Cheung and De Rosa, just proportional to the vertical (averaged to points) magnetic field

            bf = 0.5 * (bfield1 + bfield2)

            bf = 0.25*(bf[1:,1:] + bf[:-1,1:] + bf[:-1,:-1] + bf[1:,:-1])
            D = omega*(bf)

            div_test = grid.div_E(ex, ey)
            phi = ft.point_transform(-div_test + D)
            correct_x, correct_y = grid.grad(phi)
            ex += correct_x
            ey += correct_y
        
            div_test = grid.div_E(ex, ey)

            if plot:

                if np.max(np.abs(curl_test[1:-1,1:-1] + diff[1:-1,1:-1])) > 1e-10:
                    raise Exception('Electric field calculation failed')
                if np.max(np.abs(div_test[1:-1,1:-1] - D[1:-1,1:-1])) > 1e-10:
                    raise Exception('Electric field calculation failed')

                print(bfield_fname)
                print('Curl test', np.max(np.abs(curl_test[1:-1,1:-1] + diff[1:-1,1:-1])))
                print('Div Test', np.max(np.abs(div_test[1:-1,1:-1] - D[1:-1,1:-1])))
                print('____________________________________________')

            curl_test = grid.curl_E(ex, ey)

            #print('dx', grid.dx, 'dy', grid.dy)
            #Swap sign of E, as that seems to be the way forward.
            ex = -ex
            ey = -ey

            if True:   #Constrain to envelope (removes exact solution, unfortunately)

                input_xs = np.linspace(-1,1,ex.shape[0])
                input_ys = np.linspace(-1,1,ex.shape[1])

                X, Y = np.meshgrid(input_xs, input_ys, indexing = 'ij')
                edge = 0.9; steep = 20.0
                envelope = 0.5-0.5*np.tanh(np.maximum(steep*(X**2-edge), steep*(Y**2-edge)))

                envelope[-1,:] = 0.0;envelope[0,:] = 0.0;envelope[:,0] = 0.0;envelope[:,-1] = 0.0

                ex = ex*envelope

                input_xs = np.linspace(-1,1,ey.shape[0])
                input_ys = np.linspace(-1,1,ey.shape[1])

                X, Y = np.meshgrid(input_xs, input_ys, indexing = 'ij')
                edge = 0.9; steep = 20.0
                envelope = 0.5-0.5*np.tanh(np.maximum(steep*(X**2-edge), steep*(Y**2-edge)))

                envelope[-1,:] = 0.0;envelope[0,:] = 0.0;envelope[:,0] = 0.0;envelope[:,-1] = 0.0

                ey = ey*envelope

            if False:
                plt.pcolormesh(ey, cmap = 'seismic', vmin = -np.max(np.abs(ey)), vmax = np.max(np.abs(ey)))
                plt.savefig('plots/ey%d.png' % snap)
                plt.close()
                
                plt.pcolormesh(ex, cmap = 'seismic', vmin = -np.max(np.abs(ex)), vmax = np.max(np.abs(ex)))
                plt.savefig('plots/ex%d.png' % snap)
                plt.show()

            fid = netcdf_file(efield_fname, 'w')
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
            #Still doesn't seem to like it -- sums are all correct but going in with the wrong direction.
            vid = fid.createVariable('ex', 'd', ('ys','xc'))
            vid[:] = np.swapaxes(ex, 0, 1)
            vid = fid.createVariable('ey', 'd', ('yc','xs'))
            vid[:] = np.swapaxes(ey, 0, 1)
        
            fid.close()

            if initialise:  #Calculate reference helicities
                bfield_fname = '%s%04d.nc' % (data_directory, snap)
                data = netcdf_file(bfield_fname, 'r', mmap=False)
                bx = np.swapaxes(data.variables['bx'][:], 0, 1)
                by = np.swapaxes(data.variables['by'][:], 0, 1)
                bz = np.swapaxes(data.variables['bz'][:], 0, 1)
                href = compute_inplane_helicity(grid, bx, by, bz)
                hrefs.append(np.sum(href))
                trefs.append(mag_times[snap])

        if initialise:  #Calculate final one
            bfield_fname = '%s%04d.nc' % (data_directory, snap + 1)
            data = netcdf_file(bfield_fname, 'r', mmap=False)
            bx = np.swapaxes(data.variables['bx'][:], 0, 1)
            by = np.swapaxes(data.variables['by'][:], 0, 1)
            bz = np.swapaxes(data.variables['bz'][:], 0, 1)
            href = compute_inplane_helicity(grid, bx, by, bz)
            hrefs.append(np.sum(href))
            trefs.append(mag_times[snap+1])

        print('Electric fields', start, ' to ', end, ' calculated and saved.')
        if initialise:
            np.save('./hdata/h_ref.npy', np.array(hrefs))
            np.save('./hdata/t_ref.npy', np.array(trefs))
        else:
            self.halls = halls
            self.tplots = tplots
            self.omegas = omegas

