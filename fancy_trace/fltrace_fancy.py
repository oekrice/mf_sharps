#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 15:16:37 2022

@author: Oliver Rice

Python wrapper for fortran field line tracer

Uses FLH density on the photosphere to determine which lines to trace.

Needs routine from emergence_comparisons to do that. Bit of a mess...
"""


import os
import shutil
import numpy as np
import sys
from numpy import random
import time
from scipy.io import netcdf_file
import matplotlib.pyplot as plt
import pyvista as pv
import matplotlib.gridspec as gridspec
from scipy.ndimage import gaussian_filter

#from get_flh import FLH
#pv.start_xvfb()
import matplotlib
matplotlib.rcParams['text.usetex'] = True

class Fltrace():
    def __init__(self, run, snap, show = 0, allscales = []):

        if os.uname()[1] == 'brillouin.dur.ac.uk' or os.uname()[1] == 'modigliani.dur.ac.uk':
            machine_flag = 0
        elif os.uname()[1][-14:] == 'ham8.dur.ac.uk':
            machine_flag = 1
        else:
            machine_flag = 2

        machine_flag = 0

        if machine_flag == 0:
            data_directory = '/extra/tmp/trcn27/mf3d/%03d/' % run
        elif machine_flag == 1:
            data_directory = '/nobackup/trcn27/mf3d0/%03d/' % run
        elif machine_flag == 2:
            data_directory = ('../Data/%03d/' % run)
        else:
            data_directory = ('./Data/%03d/' % run)



        self.show = show
        self.machine_flag = machine_flag
        print('Data root', data_directory)
        #Establish grid parameters (can be read in from elsewhere of course)
        self.run = run
        self.snap = snap
        self.print_flag = show

        self.save_number = self.snap
        self.option = 2   #tracing a plotting options (1 for jet, 2 for emergence)

        self.data_source = 0.
        self.data_root = data_directory

        if not (self.print_flag): #In this case, assume it's fine to wait for a bit
            while not os.path.exists('%s%04d.nc' % (self.data_root, self.snap + 1)):
                time.sleep(0.1)
        else:
            pass
        #Establish random seeds for field line plotting
        if not os.path.isfile('start_seeds.txt'):
            self.start_seeds = random.rand(1000**2)
            np.savetxt('start_seeds.txt', self.start_seeds, delimiter = ',')
        else:
            self.start_seeds = np.loadtxt('start_seeds.txt', delimiter = ',')
        data = netcdf_file('%s%04d.nc' % (self.data_root, self.snap), 'r', mmap=False)
        self.bx = np.swapaxes(data.variables['bx'][:],0,2)
        self.by = np.swapaxes(data.variables['by'][:],0,2)
        self.bz = np.swapaxes(data.variables['bz'][:],0,2)
        data.close()


        #Establish start points for the field line plotting
        self.max_line_length = 100000
        self.ds = 0.05 #Tracing 'timestep' as a proportion of the grid size
        self.weakness_limit = 1e-5  #Minimum field strength to stop plotting
        self.line_plot_length = 100  #To save time while plotting, reduce the length of the plotted lines
        self.nlines = 50000#1000000#200000

        #Import bz as a test of the resolutions (and for the pyvista plot)
        self.nx = np.shape(self.bz)[0]
        self.ny = np.shape(self.bz)[1]
        self.nz = np.shape(self.bz)[2] - 1

        paras = np.loadtxt('../parameters/variables%03d.txt' % run)

        self.nx = int(paras[1])
        self.ny = int(paras[2])

        self.x0 = paras[12]; self.x1 = paras[13]
        self.y0 = paras[14]; self.y1 = paras[15]
        self.z0 = paras[16]; self.z1 = paras[17]

        self.xs = np.linspace(self.x0,self.x1,self.nx+1)
        self.ys = np.linspace(self.y0,self.y1,self.ny+1)
        self.zs = np.linspace(self.z0,self.z1,self.nz+1)

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

        self.dx = self.xs[1] - self.xs[0]
        self.dy = self.ys[1] - self.ys[0]
        self.dz = self.zs[1] - self.zs[0]

    def trace_lines(self):

        #Folder admin
        if not os.path.exists('./fl_data/'):
            os.mkdir('fl_data')

        #flh = FLH(self)    #Do field-line helicity things
        #self.flh_photo = flh.flh_photo #FLH density on the photosphere
        #self.plot_base = self.flh_photo

        self.plot_base = self.bz[:,:,0]
        #self.plot_base = np.ones((self.nx, self.ny))
        #Find start points
        self.set_starts()

        #Create runtime variables for fortran
        self.setup_tracer()
        #Do the tracing. MAY NEED TO CHANGE DATA DIRECTORY IN fltrace.f90
        self.trace_lines_fortran()
        #Plot the field lines (using pyvista)
        if not os.path.exists('./plots/'):
            os.mkdir('plots')

        #os.system('rm ./fl_data/emiss%04d.nc' % self.snap)
        #self.plot_vista()

    def plot_emiss(self,snap = -1, allscales = []):

        if snap > 0:
            self.snap = snap

        try:
            data = netcdf_file('./fl_data/emiss%04d.nc' % self.snap, 'r', mmap=False)
            print('Emissions found', self.snap)

        except:
            print('File not found')

        mag_times = np.loadtxt('../parameters/magtimes%03d.txt' % self.run)

        #self.emiss = np.swapaxes(data.variables['emiss'][:],0,2)
        self.xsum = np.swapaxes(data.variables['emiss_xsum'][:],0,1)
        self.ysum = np.swapaxes(data.variables['emiss_ysum'][:],0,1)
        self.zsum = np.swapaxes(data.variables['emiss_zsum'][:],0,1)

        xsplot = np.linspace(self.x0,self.x1,self.zsum.shape[0]+1)
        ysplot = np.linspace(self.y0,self.y1,self.zsum.shape[1]+1)
        zsplot = np.linspace(self.z0,self.z1,self.xsum.shape[1]+1)

        #fig, axs = plt.subplots(4,1, figsize = (15,7.5))
        fig = plt.figure(layout = 'constrained', figsize = (10,7.5))

        toplots = [np.sqrt(self.zsum), np.sqrt(self.ysum), np.sqrt(self.xsum)]
        ranges = [[xsplot, ysplot], [xsplot, zsplot],[ysplot, zsplot]]
        titles = ['Top', 'x face', 'y face']

        gs0 = gridspec.GridSpec(3,1,figure=  fig)


        for i in range(3):

            if i == 0:
                gs1 = gridspec.GridSpecFromSubplotSpec(1,2, subplot_spec=gs0[0])
                ax = fig.add_subplot(gs1[0])

            else:
                gs2 = gridspec.GridSpecFromSubplotSpec(1,1, subplot_spec=gs0[i])
                ax = fig.add_subplot(gs2[0])

            toplot = toplots[i]

            #Top down
            #vmax = np.percentile(toplot, 99.75)
            pc = (np.sum([toplot > 1e-5])/np.sum([toplot >= 0.0]))
            pc_total = 100 - 0.1*pc
            vmax = np.percentile(toplot, pc_total)
            vmax = max(0.1, vmax)

            if len(allscales) > 0:
                vmax = allscales[self.snap, i]
            #vmax = 1.0
            im = ax.pcolormesh(ranges[i][0], ranges[i][1],toplot.T, vmin = 0, vmax = vmax, cmap = 'inferno')
            ax.set_aspect('equal')
            ax.set_title(titles[i])


        #Then magnetogram
        ax = fig.add_subplot(gs1[1])
        toplot = self.bz[:,:,0].T
        ax.pcolormesh(self.xs, self.ys, toplot,cmap ='seismic', vmax = np.max(np.abs(toplot)), vmin = -np.max(np.abs(toplot)))
        ax.set_aspect('equal')
        ax.set_title('Lower boundary magnetogram')

        plt.suptitle('t = %03d' % mag_times[self.snap])
        plt.savefig('fancies/fancy%03d.png' % self.snap)
        if self.show:
            plt.show()
        plt.close()

    def set_starts(self):
        #Set the start positions for the lines to be traced. Will by default try to trace in both directions from this position.

        #Plot up from the surface, based on some threshold of how strong the magnetic field is... Could be fun.

        nlines = self.nlines

        nx_lines = int(np.sqrt(nlines))
        ny_lines = int(np.sqrt(nlines))

        alpha = 2.0

        self.starts = []

        for i in range(nx_lines):
             for j in range(ny_lines):
                 xstart = 0.75*self.xs[0] + (i/nx_lines)*0.75*(self.xs[-1] - self.xs[0])
                 ystart = 0.75*self.ys[0] + (j/ny_lines)*0.75*(self.ys[-1] - self.ys[0])

                 self.starts.append([xstart, ystart, 0.0])

        print('Tracing', len(self.starts), 'lines')

        self.nstarts = len(self.starts)
        self.starts = np.array(self.starts).reshape(self.nstarts*3)

    def setup_tracer(self):
        #Output runtime variables to be read-in to Fortran code
        max_line_length = 10000
        ds = 0.1 #Tracing 'timestep' as a proportion of the grid size
        print_flag = 1  #Print some things as the tracing happens
        save_all = 0

        self.nx_out = 1024   #Resolutions of the output data
        self.ny_out = int(self.nx_out*self.ny/self.nx)
        self.nz_out = int(self.nx_out*self.nz/self.nx)

        variables = np.zeros((30))

        variables[0] = self.run
        variables[1] = self.nx
        variables[2] = self.ny
        variables[3] = self.nz
        variables[4] = self.x0
        variables[5] = self.x1
        variables[6] = self.y0
        variables[7] = self.y1
        variables[8] = self.z0
        variables[9] = self.z1
        variables[10] = self.snap
        variables[11] = self.nstarts
        variables[12] = self.print_flag
        variables[13] = self.max_line_length
        variables[14] = self.ds
        variables[15] = self.weakness_limit
        variables[16] = self.data_source
        variables[17] = self.machine_flag
        variables[18] = self.nx_out
        variables[19] = self.ny_out
        variables[20] = self.nz_out
        variables[21] = print_flag
        variables[22] = save_all

        np.savetxt('./fl_data/flparameters%04d.txt' % self.snap, variables)   #variables numbered based on run number (up to 1000)
        np.savetxt('./fl_data/starts%04d.txt' % self.snap, self.starts)   #Coordinates of the start points of each field line (do this in python)

    def trace_lines_fortran(self):
        os.system('make')
        if os.uname()[1] == 'brillouin.dur.ac.uk' or os.uname()[1] == 'modigliani.dur.ac.uk':
            os.system('/usr/lib64/openmpi/bin/mpiexec -np 1 ./bin/fltrace %d' % self.snap)
        elif os.uname()[1] == 'login1.ham8.dur.ac.uk' or os.uname()[1] == 'login2.ham8.dur.ac.uk':
            os.system('mpiexec -np 1 ./bin/fltrace %d' % self.snap)
        else:
            os.system('mpirun -np 1 ./bin/fltrace %d' % self.snap)

def determine_scales(snap_min, snap_max):

    allscales = []
    for snap in range(snap_min, snap_max):

        try:
            data = netcdf_file('./fl_data/emiss%04d.nc' % snap, 'r', mmap=False)
            xsum = np.swapaxes(data.variables['emiss_xsum'][:],0,1)
            ysum = np.swapaxes(data.variables['emiss_ysum'][:],0,1)
            zsum = np.swapaxes(data.variables['emiss_zsum'][:],0,1)

            data.close()

            toplots = [np.sqrt(zsum), np.sqrt(ysum), np.sqrt(xsum)]
            scale = []
            for i, toplot in enumerate(toplots):

                pc = (np.sum([toplot > 1e-5])/np.sum([toplot >= 0.0]))
                pc_total = 100 - 0.1*pc
                vmax = np.percentile(toplot, pc_total)

                scale.append(vmax)

            allscales.append(scale)

        except:
            pass

    allscales = np.array(allscales)
    for i in range(3):
        allscales[:,i] = gaussian_filter(allscales[:,i],10)

    return allscales



if len(sys.argv) > 1:
    run = int(sys.argv[1])
else:
    run = 0

if len(sys.argv) > 2:
    snap_min = int(sys.argv[2])
else:
    snap_min = 0

if len(sys.argv) > 3:
    snap_max = int(sys.argv[3])
else:
    #How many snaps actually are there?
    paras = np.loadtxt('../parameters/variables%03d.txt' % run)

    snap_max = int(paras[5] - 1)

if len(sys.argv) > 4:
    show = int(sys.argv[4])
else:
    show = 0

print('Run', run, snap_min, snap_max)
nset = 1 #Number of concurrent runs. Receives input 0-(nset-1)

#snap_min = 0
#snap_max = 94
skip = 1

for run in range(run, run+1):
    for snap in range(snap_min, snap_max, skip):
        print('Run number', run, 'Plot number', snap_min)
        fltrace = Fltrace(run = run, snap = snap, show = show)
        fltrace.trace_lines()
        fltrace.plot_emiss()

        os.system('rm ./fl_data/flparameters%04d.txt' % fltrace.snap)
        os.system('rm ./fl_data/starts%04d.txt' % fltrace.snap)

        #snap_min = snap_min + nset

    #Redo with nicer scales...
    allscales = determine_scales(snap_min, snap_max)
    for snap in range(snap_min, min(len(allscales), snap_max), skip):
        fltrace = Fltrace(run = run, snap = snap, show = show)

        fltrace.plot_emiss(snap = snap, allscales = allscales)















