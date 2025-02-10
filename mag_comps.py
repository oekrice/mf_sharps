#New script to do helicity analysis on both the original LARE 'magnetograms' and anything generated by MF. Used to determine the 'best' value for the myserious omega.


import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib
from scipy.io import netcdf_file
from scipy.interpolate import RegularGridInterpolator
from scipy.fft import fft, ifft2, fft2, ifft
import os
from scipy.ndimage import gaussian_filter
matplotlib.rcParams['text.usetex'] = True


class Grid():
    """In the interest of doing it properly, put grid parameters in here"""
    def __init__(self, run):

        paras = np.loadtxt('parameters/variables%03d.txt' % run)

        import_resolution = 128

        #Define the grid onto which the electric field should be outputted
        self.nx = import_resolution
        self.ny = import_resolution

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

        self.xs_import = np.linspace(self.x0, self.x1, import_resolution+1)
        self.ys_import = np.linspace(self.y0, self.y1, import_resolution+1)


        self.dx_import = self.xs[1] - self.xs[0]
        self.dy_import = self.ys[1] - self.ys[0]

        self.xc_import = np.linspace(self.x0 - self.dx/2, self.x1 + self.dx/2, import_resolution+2)
        self.yc_import = np.linspace(self.y0 - self.dy/2, self.y1 + self.dy/2, import_resolution+2)


class compute_inplane_helicity():
    #Uses all three components of the magnetic field to give values for A.B at EACH of the input timesteps.
    #Requires the script from Chris' flt code
    def __init__(self, run, snap = 0):

        if snap < 0:
            #Find latest snap
            for snap_try in range(10000):
                if not os.path.exists('./mf_mags/%03d/%04d.nc' % (run, snap_try)):
                    snap = snap_try - 1
                    break

        sharp_id = 956
        grid = Grid(0)
        h_ref = []
        ts = []
        #Call run = -1 the reference case
        runs = [-1,0]
        omegas = []

        scales = np.zeros((4,2))

        for snap_num in range(snap, snap+1):
            fig, axs = plt.subplots(4, 2)

            for ri, run in enumerate(runs):
                if run >= 0:
                    paras = np.loadtxt('parameters/variables%03d.txt' % run)
                    omega = paras[28]
                    omegas.append(omega)
                    print(run, 'omega', omega)

                if (snap%10) == 0:
                    print('Run', run, ', snap', snap)
                #Find from LARE

                if run < 0:
                    source = '/extra/tmp/trcn27/sharps/%05d_mag/' % sharp_id
                else:
                    source = './mf_mags/%03d/' % run

                data_directory = source

                bfield_fname = '%s%04d.nc' % (data_directory, snap_num)

                try:
                    data = netcdf_file(bfield_fname, 'r', mmap=False)

                except:
                    print('File', bfield_fname, 'not found')
                    continue

                if run < 0:
                    bx = np.swapaxes(data.variables['bx'][:],0,1)
                    by = np.swapaxes(data.variables['by'][:],0,1)
                    bz = np.swapaxes(data.variables['bz'][:],0,1)
                    bx0 = 0.5*(bx[1:,1:-1] + bx[:-1,1:-1])
                    by0 = 0.5*(by[1:-1,1:] + by[1:-1,:-1])
                    bz0 = bz[1:-1,1:-1]

                else:
                    bx = np.swapaxes(data.variables['bx'][:],0,1)
                    by = np.swapaxes(data.variables['by'][:],0,1)
                    bz = np.swapaxes(data.variables['bz'][:],0,1)
                    bx0 = 0.5*(bx[1:,:] + bx[:-1,:])
                    by0 = 0.5*(by[:,1:] + by[:,:-1])
                    bz0 = bz[:,:]



                #Trim the edges out as these can go a bit screwy and bugger up the results

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

                #print('Vector potential test', np.max(np.abs(bz[1:-1,1:-1] - bz_test[1:-1,1:-1]))/np.max(np.abs(bz[1:-1,1:-1])))
                #This vector potential should be reasonably OK... Need code to test though

                hfield = np.sqrt(np.abs(ax0*bx0 + ay0*by0 + az0*bz0))

                print(np.sum(hfield))
                toplots = [bx.T, by.T, bz.T, hfield.T]

                for plot in range(4):
                    if ri == 0:
                        vmin = -np.max(np.abs(toplots[plot]))
                        vmax = np.max(np.abs(toplots[plot]))
                        scales[plot, 0] = vmin
                        scales[plot, 1] = vmax

                    else:
                        vmin = scales[plot,0]
                        vmax = scales[plot,1]
                    ax = axs[plot,ri]
                    im = ax.pcolormesh(toplots[plot], cmap = 'seismic', vmin = vmin, vmax = vmax)

                    plt.colorbar(im, ax = ax)

                if ri > 0:
                    plt.tight_layout()
                    plt.show()

                if False:
                    string = ['Reference', 'MF'][ri]
                    if run < 0:
                        fig, axs = plt.subplots(2,4, figsize = (10,5))

                    if run < 0:
                        row = 0

                    else:
                        row = ri - 1

                    ax = axs[ri, 0]
                    im = ax.pcolormesh(bx.T, cmap = 'seismic', vmax = np.max(np.abs(bx)), vmin = -np.max(np.abs(bx)))
                    plt.colorbar(im, ax = ax)
                    ax.set_title('$B_x$, %s' % string)

                    ax = axs[ri, 1]
                    im =axs[ri,1].pcolormesh(by.T, cmap = 'seismic', vmax = np.max(np.abs(by)), vmin = -np.max(np.abs(by)))
                    plt.colorbar(im, ax = ax)
                    ax.set_title('$B_y$, %s' % string)

                    ax = axs[ri, 2]
                    im = axs[ri,2].pcolormesh(bz.T, cmap = 'seismic', vmax = np.max(np.abs(bz)), vmin = -np.max(np.abs(bz)))
                    plt.colorbar(im, ax = ax)
                    ax.set_title('$B_z$, %s' % string)

                    ax = axs[ri, 3]
                    im = axs[ri,3].pcolormesh(hfield.T, cmap = 'seismic', vmax = np.max(np.abs(hfield)), vmin = -np.max(np.abs(hfield)))
                    plt.colorbar(im, ax = ax)
                    ax.set_title('$(A.B)$, %s' % string)

                hfield = np.sqrt(np.abs(hfield))

            plt.suptitle('t = %d' % (snap*0.5))
            plt.tight_layout()
            #plt.savefig('./hplots/%04d.png' % snap)
            #plt.show()
            plt.close()

        #np.save('./hdata/h_ref.npy', h_ref)
        #np.save('./hdata/ts.npy', ts)

        if False:
            fig = plt.figure(figsize = (10,7))
            for ri in range(len(runs)-1):
                plt.plot(ts[:len(h_all[ri])], h_all[ri], label = ('Omega factor= %.4f' % omegas[ri]))

            plt.plot(ts, h_ref, c= 'black', linestyle = 'dashed', label = 'LARE Reference')

            plt.legend()

            plt.xlabel('Time')
            plt.ylabel('Surface helicity A.B')
            plt.tight_layout()
            plt.savefig('helicityanalysis.png')
            plt.show()



if len(sys.argv) > 1:
    snap = int(sys.argv[1])
else:
    snap = -1

compute_inplane_helicity(0, snap = snap)


