import numpy as np
import matplotlib.pyplot as plt
import random
from scipy.linalg import eigh_tridiagonal
from scipy.io import netcdf_file, FortranFile
from scipy.interpolate import RegularGridInterpolator


class compute_initial_condition():
    def __init__(self, grid, lbound_fn, run, background_strength = 0.0, background_angle = 0.0, boundary_error_limit = 0.0, init_filename = '/inits/init000.nc'):
        #use the same grid notation as in the proper diagnostic calculator, even though it can be bodged for now

        print('Calculating initial potential field for given boundary function...')
        self.xs = np.linspace(grid.x0,grid.x1,grid.nx+1)
        self.ys = np.linspace(grid.y0,grid.y1,grid.ny+1)
        self.zs = np.linspace(grid.z0,grid.z1,grid.nz+1)

        self.nx = grid.nx
        self.ny = grid.ny
        self.nz = grid.nz

        self.xc = np.zeros(self.nx + 2)
        self.yc = np.zeros(self.ny + 2)
        self.zc = np.zeros(self.nz + 2)

        self.xc[1:-1] = 0.5*(self.xs[1:] + self.xs[:-1])
        self.yc[1:-1] = 0.5*(self.ys[1:] + self.ys[:-1])
        self.zc[1:-1] = 0.5*(self.zs[1:] + self.zs[:-1])

        self.xc[0] = self.xc[1] - (self.xc[2] - self.xc[1])
        self.yc[0] = self.yc[0] - (self.yc[2] - self.yc[2])
        self.zc[0] = self.zc[0] - (self.zc[2] - self.zc[2])

        self.xc[-1] = self.xc[-2] + (self.xc[-2] - self.xc[-3])
        self.yc[-1] = self.yc[-2] + (self.yc[-2] - self.yc[-3])
        self.zc[-1] = self.zc[-2] + (self.zc[-2] - self.zc[-3])


        self.dx = np.sum(self.xs[1:] - self.xs[:-1])/len(self.xs[1:])
        self.dy = np.sum(self.ys[1:] - self.ys[:-1])/len(self.ys[1:])
        self.dz = np.sum(self.zs[1:] - self.zs[:-1])/len(self.zs[1:])

        self.init_filename = init_filename
        self.boundary_error_limit = boundary_error_limit
        print('Initialising with potential field')

        #Need to interpolate onto a new grid. OR NOT
        self.lbound = lbound_fn[1:-1,1:-1]

        self.lbound_transform = np.zeros((self.nx,self.ny))

        self.xbasis = np.zeros((self.nx, self.nx+2))
        self.ybasis = np.zeros((self.ny, self.ny+2))

        #find eigenvalues and basis vectors (eigenvalues are m^2 and numbered k)
        #-----------------------------------------
        print('Calculating horizontal basis functions...')
        self.m2, self.n2, self.xbasis[:,1:-1], self.ybasis[:,1:-1] = self.find_eigenstuff()

        self.xbasis[:,0] = self.xbasis[:,1]; self.xbasis[:,-1] = self.xbasis[:,-2]
        self.ybasis[:,0] = self.ybasis[:,1]; self.ybasis[:,-1] = self.ybasis[:,-2]

        print('Starting "Fourier" Transform')
        error = 1e6
        self.phi = np.zeros((self.nx+2, self.ny+2,self.nz+2))
        self.test1 = np.zeros((self.nx+2))
        for k_c in range(self.nx+self.ny):
            for k1 in range(self.nx):
                k2 = k_c-k1
                if k2 >= self.ny or k2 < 0:
                    continue
                if error < boundary_error_limit:
                    break

                self.lbound_transform[k1,k2] = self.coeff(self.lbound, k1,k2)
                if abs(self.lbound_transform[k1,k2]) < 1e-10:
                    continue
                zbasis = self.find_zbasis(self.m2[k1],self.n2[k2])
                self.phi = self.phi + self.lbound_transform[k1,k2]*self.xbasis[k1,:][:,np.newaxis,np.newaxis]*self.ybasis[k2,:][np.newaxis,:,np.newaxis]*zbasis[np.newaxis,np.newaxis,:]
                self.test1 = self.test1 + self.lbound_transform[k1,k2]*self.xbasis[k1,:]*self.ybasis[k2,0]
                lbound_test = (self.phi[1:-1,1:-1,1] - self.phi[1:-1,1:-1,0])/self.dz
                error = np.sqrt(np.sum((lbound_test - self.lbound)**2)/(self.nx*self.ny))

        self.bx = np.zeros((self.nx+1,self.ny+2,self.nz+2))
        self.by = np.zeros((self.nx+2,self.ny+1,self.nz+2))
        self.bz = np.zeros((self.nx+2,self.ny+2,self.nz+1))

        self.bx[:,1:-1,1:-1] = (self.phi[1:,1:-1,1:-1] - self.phi[:-1,1:-1,1:-1])/self.dx
        self.by[1:-1,:,1:-1] = (self.phi[1:-1,1:,1:-1] - self.phi[1:-1,:-1,1:-1])/self.dy
        self.bz[1:-1,1:-1,:] = (self.phi[1:-1,1:-1,1:] - self.phi[1:-1,1:-1,:-1])/self.dz

        #Jx Conditions
        self.by[:,:,0] = self.by[:,:,1] - self.dz*(self.bz[:,1:,0] - self.bz[:,:-1,0])/self.dy
        self.by[:,:,-1] = self.by[:,:,-2] + self.dz*(self.bz[:,1:,-1] - self.bz[:,:-1,-1])/self.dy

        self.bz[:,0,:] = self.bz[:,1,:] - self.dy*(self.by[:,0,1:] - self.by[:,0,:-1])/self.dz
        self.bz[:,-1,:] = self.bz[:,-2,:] + self.dy*(self.by[:,-1,1:] - self.by[:,-1,:-1])/self.dz

        #Jy Conditions
        self.bx[:,:,0] = self.bx[:,:,1] - self.dz*(self.bz[1:,:,0] - self.bz[:-1,:,0])/self.dx
        self.bx[:,:,-1] = self.bx[:,:,-2] + self.dz*(self.bz[1:,:,-1] - self.bz[:-1,:,-1])/self.dx

        self.bz[0,:,:] = self.bz[1,:,:] - self.dx*(self.bx[0,:,1:] - self.bx[0,:,:-1])/self.dz
        self.bz[-1,:,:] = self.bz[-2,:,:] + self.dx*(self.bx[-1,:,1:] - self.bx[-1,:,:-1])/self.dz

        #Jz Conditions
        self.by[0,:,:] = self.by[1,:,:] - self.dx*(self.bx[0,1:,:] - self.bx[0,:-1,:])/self.dy
        self.by[-1,:,:] = self.by[-2,:,:] + self.dx*(self.bx[-1,1:,:] - self.bx[-1,:-1,:])/self.dy

        self.bx[:,0,:] = self.bx[:,1,:] - self.dy*(self.by[1:,0,:] - self.by[:-1,0,:])/self.dx
        self.bx[:,-1,:] = self.bx[:,-2,:] + self.dy*(self.by[1:,-1,:] - self.by[:-1,-1,:])/self.dx

        print('')
        print('Field Calculated')

        self.test_phi()

        self.find_vector_potentials()

        self.save_initial_condition()


    def find_eigenstuff(self):
        # Uses scipy tridiagonal solver to find the numerical approximations to the sine functions that have the desired properties.
        #Generates a matrix etc. then solves. Should only need to do this once for a given resolution. Doesn't depend on boundary conditions etc.

        #x FIRST
        d = 2*np.ones(self.nx)
        e = -1*np.ones(self.nx-1)
        d[0] = 1.0; d[-1] = 1.0
        w, xb = eigh_tridiagonal(d, e)
        m2 = w/(self.dx**2)
        #THEN y
        d = 2*np.ones(self.ny)
        e = -1*np.ones(self.ny-1)
        d[0] = 1.0; d[-1] = 1.0
        w, yb = eigh_tridiagonal(d, e)
        n2 = w/(self.dy**2)

        return m2, n2, xb.T, yb.T

    def find_zbasis(self, m2, n2):
        zbasis = np.zeros((self.nz+2))  #number of x modes plus dimension (with ghosts) in the y direction
        zbasis[-1] = 1.0   #set to zero at the top (roughly)
        zbasis[-2] = -1.0
        for i in range(self.nz, 0, -1):
            zbasis[i-1] = (m2+n2)*zbasis[i]*self.dz**2
            zbasis[i-1] += 2*zbasis[i] - zbasis[i+1]
        dfact = (zbasis[1] - zbasis[0])/self.dz
        zbasis = zbasis/dfact
        return zbasis

    def coeff(self, lbound, k1, k2):

        #Finds the coefficient of the k1, k2 mode
        return np.sum(lbound[:,:]*self.xbasis[k1,1:-1][:,np.newaxis]*self.ybasis[k2,1:-1][np.newaxis,:])/np.sum(self.xbasis[k1,1:-1][:,np.newaxis]**2*self.ybasis[k2,1:-1][np.newaxis,:]**2)

    def find_az(self):
        #Finds the vector potential from the magnetic field
        print(np.shape(self.bxp), np.shape(self.byp))
        az = np.zeros((self.nx+1,self.ny+1))
        for i in range(1,self.nx+1):
            az[i,0] = az[i-1,0] + self.dx*self.byp[i-1,0]
            for j in range(1,self.ny+1):
                az[i,j] = az[i,j-1] - self.dy*self.bxp[i,j-1]
        return az

    def find_vector_potentials(self):
        #Finds the vector potentials ax, ay, az from the existing magnetic field
        #To be ported to fortran to ensure each core everything matches up nicely.
        #Bottom left corner set as zero, for now at least
        self.ax = np.zeros((self.nx, self.ny+1, self.nz+1))
        self.ay = np.zeros((self.nx+1, self.ny, self.nz+1))
        self.az = np.zeros((self.nx+1, self.ny+1, self.nz))

        #Do Ax first as one would do in 2D, just on the bottom
        for i in range(0,self.nx):
            for j in range(1,self.ny+1):
                self.ax[i,j,0] = -self.bz[i+1,j,0]*self.dy + self.ax[i,j-1,0]
        #Go up through the layers
        for k in range(1,self.nz+1):
            #Ax first
            for i in range(0,self.nx):
                for j in range(0,self.ny+1):
                    self.ax[i,j,k] = self.ax[i,j,k-1] + self.dz*self.by[i+1,j,k]
            #Then Ay
            for i in range(0,self.nx+1):
                for j in range(0,self.ny):
                    self.ay[i,j,k] = self.ay[i,j,k-1] - self.dz*self.bx[i,j+1,k]

        #Add random fluctuations
        self.ax = self.ax + 0.0*np.random.random(size = np.shape(self.ax))
        self.ay = self.ay + 0.0*np.random.random(size = np.shape(self.ay))
        self.az = self.az + 0.0*np.random.random(size = np.shape(self.az))

        bx_test = np.zeros((self.nx+1,self.ny+2,self.nz+2))
        by_test = np.zeros((self.nx+2,self.ny+1,self.nz+2))
        bz_test = np.zeros((self.nx+2,self.ny+2,self.nz+1))

        bx_test[:,1:-1,1:-1] = (self.az[:,1:,:] - self.az[:,:-1,:])/self.dy - (self.ay[:,:,1:] - self.ay[:,:,:-1])/self.dz
        by_test[1:-1,:,1:-1] = (self.ax[:,:,1:] - self.ax[:,:,:-1])/self.dz - (self.az[1:,:,:] - self.az[:-1,:,:])/self.dx
        bz_test[1:-1,1:-1,:] = (self.ay[1:,:,:] - self.ay[:-1,:,:])/self.dx - (self.ax[:,1:,:] - self.ax[:,:-1,:])/self.dy

        if False:
            print('Vector potential errors:')
            print(np.max(np.abs(bx_test[:,1:-1,1:-1] - self.bx[:,1:-1,1:-1])))
            print(np.max(np.abs(by_test[1:-1,:,1:-1] - self.by[1:-1,:,1:-1])))
            print(np.max(np.abs(bz_test[1:-1,1:-1,:] - self.bz[1:-1,1:-1,:])))

            print('Lbound Flux', np.sum(np.abs(bz_test[:,:,0]))*(self.dx*self.dy))
            print('Open Flux', np.sum(np.abs(bz_test[:,:,-1]))*(self.dx*self.dy))
            print('----------------------------------------------')

    def second_derivative(self, array, d):
        return (array[:-2] - 2*array[1:-1] + array[2:])/d**2

    def fcheck(self, bound_trans):
        bcheck = 0.0*ubound
        for k in range(self.nx -1):
            bcheck = bcheck + bound_trans[k]*self.xbasis[k,1:-1]
        return bcheck

    def mode(self, m):
        return np.sin(0.5*m*np.pi*self.xc/self.xs[-1])

    def test_phi(self):
        d2x = (self.phi[:-2,1:-1,1:-1] - 2*self.phi[1:-1,1:-1,1:-1] + self.phi[2:,1:-1,1:-1])/self.dx**2
        d2y = (self.phi[1:-1,:-2,1:-1] - 2*self.phi[1:-1,1:-1,1:-1] + self.phi[1:-1,2:,1:-1])/self.dy**2
        d2z = (self.phi[1:-1,1:-1,:-2] - 2*self.phi[1:-1,1:-1,1:-1] + self.phi[1:-1,1:-1,2:])/self.dz**2
        print('----------------------------------------------')
        print('Max laplacian', np.max(np.abs(d2x + d2y + d2z)))
        lbound_test = (self.phi[1:-1,1:-1,1] - self.phi[1:-1,1:-1,0])/self.dz
        print('Max lbound error', np.max(np.abs(self.lbound - lbound_test)))

        jx =  (self.bz[2:-2,2:-1,1:-1] - self.bz[2:-2,1:-2,1:-1])/self.dy - (self.by[2:-2,1:-1,2:-1] - self.by[2:-2,1:-1,1:-2])/self.dz
        jy =  (self.bx[1:-1,2:-2,2:-1] - self.bx[1:-1,2:-2,1:-2])/self.dz - (self.bz[2:-1,2:-2,1:-1] - self.bz[1:-2,2:-2,1:-1])/self.dx
        jz =  (self.by[2:-1,1:-1,2:-2] - self.by[1:-2,1:-1,2:-2])/self.dx - (self.bx[1:-1,2:-1,2:-2] - self.bx[1:-1,1:-2,2:-2])/self.dy
        print('Max interior currents', np.max(np.abs(jx)),np.max(np.abs(jy)),np.max(np.abs(jz)))
        print('----------------------------------------------')

        #self.bx = np.zeros((self.nx+1,self.ny+2,self.nz+2))
        #self.by = np.zeros((self.nx+2,self.ny+1,self.nz+2))
        #self.bz = np.zeros((self.nx+2,self.ny+2,self.nz+1))

        #Test currents (including those on the boundaries, as required for A update)
        jx = (self.bz[1:-1,1:,:] - self.bz[1:-1,:-1,:])/self.dy - (self.by[1:-1,:,1:] - self.by[1:-1,:,:-1])/self.dz
        jy =  (self.bx[:,1:-1,1:] - self.bx[:,1:-1,:-1])/self.dz - (self.bz[1:,1:-1,:] - self.bz[:-1,1:-1,:])/self.dx
        jz =  (self.by[1:,:,1:-1] - self.by[:-1,:,1:-1])/self.dx - (self.bx[:,1:,1:-1] - self.bx[:,:-1,1:-1])/self.dy

        print('Max currents', np.max(np.abs(jx)),np.max(np.abs(jy)),np.max(np.abs(jz)))
        print('----------------------------------------------')

        #If errors are too large, abort
        if np.sqrt(np.sum((lbound_test - self.lbound)**2)/(self.nx*self.ny)) > self.boundary_error_limit:
            raise Exception('Potential field calculation failed')
        if np.max(np.abs(jx)) + np.max(np.abs(jy)) + np.max(np.abs(jz)) > 1e-6:
            raise Exception('Potential field calculation failed')


    def save_initial_condition(self):
        #Saves the magnetic field to the netcdf file in the inits folder
        #Can also save the vector potential (in some gauge) if I'm smart enough to do that
        #Save data out as netcdf file. Taking care with dimensions
        fid = netcdf_file(self.init_filename, 'w')
        fid.createDimension('xs', self.nx+1)
        fid.createDimension('xc', self.nx)
        fid.createDimension('ys', self.ny+1)
        fid.createDimension('yc', self.ny)
        fid.createDimension('zs', self.nz+1)
        fid.createDimension('zc', self.nz)

        vid = fid.createVariable('xs', 'd', ('xs',))
        vid[:] = self.xs
        vid = fid.createVariable('xc', 'd', ('xc',))
        vid[:] = self.xc[1:-1]
        vid = fid.createVariable('ys', 'd', ('ys',))
        vid[:] = self.ys
        vid = fid.createVariable('yc', 'd', ('yc',))
        vid[:] = self.yc[1:-1]
        vid = fid.createVariable('zs', 'd', ('zs',))
        vid[:] = self.zs
        vid = fid.createVariable('zc', 'd', ('zc',))
        vid[:] = self.zc[1:-1]

        #Transposes are necessary as it's easier to flip here than in Fortran
        vid = fid.createVariable('bx', 'd', ('zc','yc','xs'))
        vid[:] = self.bx[:,1:-1,1:-1].T
        vid = fid.createVariable('by', 'd', ('zc','ys','xc'))
        vid[:] = self.by[1:-1,:,1:-1].T
        vid = fid.createVariable('bz', 'd', ('zs','yc','xc'))
        vid[:] = self.bz[1:-1,1:-1,:].T

        vid = fid.createVariable('ax', 'd', ('zs','ys','xc'))
        vid[:] = self.ax.T
        vid = fid.createVariable('ay', 'd', ('zs','yc','xs'))
        vid[:] = self.ay.T
        vid = fid.createVariable('az', 'd', ('zc','ys','xs'))
        vid[:] = self.az.T

        fid.close()
































