#This script is designed such that the magnetofriction can ITERATIVELY pick the best omega to match the lower boundary nicely at all times.

#Also designed to integrate downloading the SHARPS data for an active region and figuring all that stuff out.

#SHOULD do everything automatically. This could be fun...
#Integrate with fancy_trace too?

import os
import shutil
import numpy as np
import sys
from numpy import random
import time
import matplotlib.pyplot as plt

from init import compute_initial_condition
from write_electric import compute_electrics, compute_inplane_helicity, read_boundary
from scipy.io import netcdf_file
from scipy.optimize import curve_fit
from scipy.io import netcdf_file

from obtain_sharps import obtain_sharp, lookup_sharp
from convert_sharps import sharp_info, convert_sharp



if len(sys.argv) > 1:
    run = int(sys.argv[1])
else:
    run = 0

if len(sys.argv) > 2:
    nprocs = int(sys.argv[2])
else:
    nprocs = 1

if len(sys.argv) > 3:
    region_id = int(sys.argv[3])
else:
    region_id = 11318

try:
    if os.uname()[1] == 'brillouin.dur.ac.uk':
        hflag = 0
        winflag = 0
    else:
        hflag = 1
        winflag = 0
except:
    hflag = 0
    winflag = 1

sharp_id = 956 #1449

sharps_directory = '/extra/tmp/trcn27/sharps/'
max_mags = 100 #Maximum number of input magnetograms (won't convert all the import data)
time_per_snap = 0.05  #Time units per input minute

mag_start = 0   #First magnetogram to start from

normalise_inputs = True       #If True, will normalise all the magnetic fields such that the max radial component is 1. Also adresses flux balance.
dothings = False
check_data = False
recalculate_inputs = dothings   #Redo the interpolation from the SHARP inputs onto this grid
recalculate_init = dothings       #Recalculates the initial potential field
recalculate_boundary = True  #Recalculates the initial boundary conditions (zero-Omega) and the reference helicity

nx = 128

#DYNAMIC SYSTEM PARAMETERS
#-------------------------------------
voutfact = 0.0   #Outflow speed. If negative, will go for as much as possible without instabilities
shearfact = 0.0#3.7e-5   #factor by which to change the imported 'speed'
eta0 = 0.0

tstart = 0.0

ndiags = 1000
nplots = -1

nu0 = 0.0
eta = 1.0#5e-4*nu0
#eta = 10.0

x0 = -100.0; x1 = 100.0   #Keep this as 100, no matter what the domain size is. Because that seemed to work.
y0 = -100.0; y1 = 100.0
z0 = 0.0; z1 = 100.0

init_number = run
omega = 0.0

#Variables for the pressure term
decay_type = 0  #Decay types -- 0 for none, 1 for exponential, 2/3 for tanh. Same as the 2D cases.

if decay_type == 0: #No pressure
    zstar = 0.0; a = 0.0; b = 0.0; deltaz = 0.0

if decay_type == 1: #exponential decay
    zstar = run#np.linspace(0.0,0.3,10)[run//50]*z1
    if zstar > 0:
        b = zstar/np.log(2)
        a = 0.5
        deltaz = 0.0*z1
    else:
        decay_type = 0
        zstar = 0.0; a = 0.0; b = 0.0; deltaz = 0.0

if decay_type == 2: #smooth tanh
    a = 0.25; b = 1.0
    zstar = np.linspace(0.0,0.3,10)[run]*z1
    deltaz = 0.1*z1

if decay_type == 3: #sharp tanh
    a = 0.25; b = 1.0
    zstar = np.linspace(0.0,0.3,7)[run]*z1
    deltaz = 0.02*z1

#SOME FOLDER ADMIN
#-------------------------------------

if not hflag:
    data_directory = '/extra/tmp/trcn27/mf3d/%03d/' % run
else:
    data_directory = '/nobackup/trcn27/mf3d0/%03d/' % run

if winflag:
    data_directory = '/Data/'
    
if not os.path.isdir('./inits'):
    os.mkdir('./inits')

if not os.path.isdir('./parameters'):
    os.mkdir('./parameters')

if not os.path.isdir('./diagnostics'):
    os.mkdir('./diagnostics')

if not os.path.isdir('./hdata'):
    os.mkdir('./hdata')

if not os.path.isdir(data_directory[:-4]):
    os.mkdir(data_directory[:-4])

if os.path.isdir(data_directory) and mag_start == 0:
    for i in range(1000):
        if os.path.isfile('%s%04d.nc' % (data_directory, i)):
            os.remove('%s%04d.nc' % (data_directory, i))
elif not os.path.isdir(data_directory):
    os.mkdir(data_directory)

if not os.path.isdir('./mf_mags/'):
    os.mkdir('./mf_mags/')

#Make magnetogram filenames if necessary
if os.path.isdir('./mf_mags/%03d/' % run) and mag_start == 0:
    for i in range(1000):
        if os.path.isfile('./mf_mags/%03d/%04d.nc' % (run, i)):
            os.remove('./mf_mags/%03d/%04d.nc' % (run, i))
elif not os.path.isdir('./mf_mags/%03d/' % run):
    os.mkdir('./mf_mags/%03d/' % run)

#Create initial condition using new init.py (potential field with arbitrary lower boundary and domain dimensions)
#-------------------------------------

class Grid():
    def __init__(self):
        self.x0 = x0; self.x1 = x1
        self.y0 = y0; self.y1 = y1
        self.z0 = z0; self.z1 = z1
        self.nx = nx ; self.ny = ny; self.nz = nz
        self.dx = (self.x1 - self.x0)/nx
        self.xs = np.linspace(self.x0,self.x1,self.nx+1)
        self.xc = np.linspace(self.x0-self.dx/2, self.x1+self.dx/2, self.nx+2)

        #Alter these ranges as appropriate based on the new resolutions
        self.y0 = self.x0*self.ny/self.nx
        self.y1 = self.x1*self.ny/self.nx
        self.dy = (self.y1 - self.y0)/ny
        self.ys = np.linspace(self.y0,self.y1,self.ny+1)
        self.yc = np.linspace(self.y0-self.dy/2, self.y1+self.dy/2, self.ny+2)

        self.z0 = self.z0*self.nz/self.nx
        self.z1 = self.z1*self.nz/self.nx
        self.dz = (self.z1 - self.z0)/nz
        self.zs = np.linspace(self.z0,self.z1,self.nz+1)
        self.zc = np.linspace(self.z0-self.dz/2, self.z1+self.dz/2, self.nz+2)


#Step 0: If sharp data isn't downloaded yet, do that. May take some time.

if sharp_id < 0:
    sharp_id = lookup_sharp(region_id)   #Attempt to find valid HARP number from the active region id
    print('SHARP ID found:', sharp_id)

if check_data or not os.path.exists(sharps_directory + '%05d_raw/' % sharp_id):
    if not os.path.exists(sharps_directory + '%05d_raw/' % sharp_id):
        print('No data exists. Will attempt to download...')
        obtain_sharp(sharp_id, sharps_directory, plot_flag = True)
    else:
        print('Some raw sharp data found. Checking if this is everything...')
        if os.path.exists('./parameters/raw_times%05d.npy' % sharp_id):
            #Check number of mags from the raw times
            nmags = len(np.load('./parameters/raw_times%05d.npy' % sharp_id))
            print(nmags)
            if os.path.exists(sharps_directory + '%05d_raw/%05d.nc' % (sharp_id, nmags-1)):
                print('Data exists up to the end of the required time')
            else:
                print('Data download incomplete. Trying to download remaining files...')
                obtain_sharp(sharp_id, sharps_directory, plot_flag = True)
        else:
            #Doesn't exist yet, so need to (slowly) obtain the data for how many mags
            print('Data somehow incomplete. Trying to download again...')
            obtain_sharp(sharp_id, sharps_directory, plot_flag = True)


print('____________________________________________')
#Step 1: Import sharp magnetograms and do some processing on them.
#_____________________________________________________________________

#Takes the raw .nc files and reduces (using 2D smoothing) to a sensible resolution
#Will save these out, for now, but should be deleted as per.
#Also onto staggered grid so they match with the original synthesised magnetograms
print('Importing SHARP data from number', sharp_id)

import_nx, import_ny, start, end = sharp_info(sharp_id, sharps_directory)

print('Using raw data from ', start, ' to ', end)
#Use the import ratios to establish the grid.
#Resolutions will be based on the ratio of the imported data, based on the nx above
aspect = import_ny/import_nx
ny = int(nx*aspect)
ny = 4*(ny//4)
nz = min(nx, ny)

print('Variables established. Making grid.')
print('New grid resolutions:', nx, ny, nz)

print('____________________________________________')

grid = Grid()

#Check if converted SHARPs already exist, to avoid doing this unecessarily
if os.path.exists(sharps_directory + '%05d_mag/' % sharp_id):
    if len(os.listdir(sharps_directory + '%05d_mag/' % sharp_id)) > 1 and not recalculate_inputs:
        print('Importing existing converted magnetic field...')
    else:
        print('Converting raw SHARPS onto this grid')
        convert_sharp(grid, sharp_id, sharps_directory, start=start, end=end, max_mags = max_mags, plot = False, normalise = normalise_inputs)

else:
    print('Converting raw SHARPS onto this grid')
    convert_sharp(grid, sharp_id, sharps_directory, start=start, end=end, max_mags = max_mags, plot = False, normalise = normalise_inputs)

nmags = len(os.listdir(sharps_directory + '%05d_mag/' % sharp_id))

if os.path.exists('./parameters/raw_times%05d.npy' % sharp_id):
    raw_mag_times = np.load('./parameters/raw_times%05d.npy' % sharp_id)
else:
    raw_mag_times = np.linspace(0.0, last_magnetogram_time, nmags)

mag_times = np.load('./parameters/mag_times%05d.npy' % sharp_id)
mag_times = (mag_times - mag_times[0])*time_per_snap
tmax = mag_times[-1]

nplots = nmags

print(nmags, 'input magnetograms found')
print('____________________________________________')

#SAVE VARIABLES TO BE READ BY FORTRAN. THESE SHOULD BE CONSTANT ON A SINGLE RUN
#-------------------------------------
#MAG MIN AND MAX CAN BE CHANGED FOR THE INDIVIDUAL RUNS LATER

variables = np.zeros((40))
#Basics
variables[0] = run
variables[1] = nx
variables[2] = ny
variables[3] = nz
variables[4] = tmax
#Outputs
variables[5] = nplots
variables[6] = ndiags
#Parameters
variables[7] = voutfact
variables[8] = shearfact
variables[9] = eta
variables[10] = nu0
variables[11] = eta0
#Grid things
variables[12] = grid.x0
variables[13] = grid.x1
variables[14] = grid.y0
variables[15] = grid.y1
variables[16] = grid.z0
variables[17] = grid.z1

#Pressure things
variables[18] = a
variables[19] = b
variables[20] = deltaz
variables[21] = zstar

variables[22] = hflag

variables[23] = decay_type

variables[24] = 0.0

#Number of imported magnetograms
variables[25] = nmags
variables[26] = tstart

variables[27] = init_number #Code to give the magnetograms and electric fields so they don't need to be done every time'

variables[28] = omega   #Just for plotting and things. Not used in the Fortran

np.savetxt('parameters/variables%03d.txt' % run, variables)   #variables numbered based on run number (up to 1000)
np.savetxt('parameters/magtimes%03d.txt' % run, mag_times)   #variables numbered based on run number (up to 1000)

#_______________________________________________________________________________________

#All the above needs to be set up once
print('Initial mag to build field from = ', mag_start)

#Find lbound from magnetogram file
mag_fname = sharps_directory + '%05d_mag/%04d.nc' % (sharp_id, mag_start)

try:
    data = netcdf_file(mag_fname, 'r', mmap=False)
    print('File', mag_fname, 'found')

except:
    print('File', mag_fname, 'not found')
    raise Exception('Initial boundary data not found')

bz = np.swapaxes(data.variables['bz'][:],0,1)

#Compute initial condition at the start -- but not later on as this will be obtained from a snapshot

if recalculate_init or not os.path.exists('./inits/init%03d.nc' % init_number):
    print('Calculating initial condition...')
    init = compute_initial_condition(grid, bz, run, boundary_error_limit = 1e-6, init_filename = './inits/init%03d.nc' % init_number)
else:
    #Put a check in to make sure the resolutions are correct
    data = netcdf_file('./inits/init%03d.nc' % init_number, 'r', mmap=False)
    init_res = np.swapaxes(data.variables['bz'][:], 0, 2).shape
    data.close()

    if init_res[0] == grid.nx and init_res[1] == grid.ny and init_res[2] == grid.nz + 1:
        print('Using existing initial condition, which appears to be the correct resolution.')
    else:
        raise Exception('Initial condition is the wrong resolution -- need to recalculate')

print('____________________________________________')

print('Calculating initial boundary conditions...')

mag_root = sharps_directory + '%05d_mag/' % sharp_id

#Compute lower boundary electrics for a zero omega, initially.
#This will save out reference helicities, as well.
if recalculate_boundary:
    compute_electrics(run, init_number, mag_root, mag_times, omega = 0.0, start = 0, end = nmags-1, initialise = True, plot = False)

stop

bx, by, bz = read_boundary('./inits/init%03d.nc' % init_number)
check = np.sum(compute_inplane_helicity(grid, bx, by, bz))

if mag_start != 0:
    halls = np.load('./hdata/halls%03d.npy' % run)
    omegas = np.load('./hdata/omegas%03d.npy' % omegas)
    tplots = np.load('./hdata/tplots%03d.npy' % tplots)
    halls[mag_start+1:] = 0.0
    omegas[mag_start+1:] = 0.0
    tplots[mag_start+1:] = 0.0

else:
    halls = np.zeros(nmags)
    omegas = np.zeros(nmags)
    tplots = np.zeros(nmags)
    halls[0] = check
    omegas[0] = omega
    tplots[0] = mag_times[0]

if hflag < 0.5:
    os.system('make')

hrefs = np.load('./hdata/h_ref.npy').tolist()

nmags_per_run = 1   #How many magnetic field INTERVALS to run for a given magnetofrictional chunk, with constant omega in each case

for block_start in range(mag_start, nmags-1, nmags_per_run):#nmags-1, nmags_per_run):

    if block_start > 0:
        omega = omegas[block_start-1]
    block_end = min(block_start + nmags_per_run, nmags-1)

    xs = []; ys = []  #For the function interpolation

    go  = True
    stopnext = False

    hardmax = 0.1
    while go:
        go = False
        omega = 0.0
        #Find the ideal omega
        if stopnext:
            go = False

        if mag_start > 0:
            back = max(0, block_start - 5)
            maxomega = max(np.abs(omegas[block_start-back:block_start])*2.0)  #Change this dynamically to keep things within range
        else:
            maxomega = hardmax

        maxomega = min(maxomega, hardmax)
        minomega = -maxomega

        omega = max(minomega, omega)
        omega = min(maxomega, omega)

        omega_range = maxomega - minomega
        print('Running step from', block_start, 'to', block_end, 'omega = ', omega)
        print('Minmax', minomega, maxomega)
        efield_data = compute_electrics(run, init_number, mag_root, mag_times, omega = omega, start = block_start, end = block_end, initialise = False, plot = False)

        variables[29] = block_start
        variables[30] = block_end

        np.savetxt('parameters/variables%03d.txt' % run, variables)   #variables numbered based on run number (up to 1000)

        #print('Using output directory "%s"' % (data_directory))
        if nprocs <= 4:
            os.system('/usr/lib64/openmpi/bin/mpiexec ffpe-summary=none -np %d ./bin/mf3d %d' % (nprocs, run))
        else:
            os.system('/usr/lib64/openmpi/bin/mpiexec -np %d --oversubscribe ./bin/mf3d %d' % (nprocs, run))

        #Check helicity against reference

        #Want to add the option to use other measures here, so probably shouldn't call it helicity
        end_fname = './mf_mags/%03d/%04d.nc' % (run, block_end)

        bx, by, bz = read_boundary(end_fname)
        check = np.sum(compute_inplane_helicity(grid, bx, by, bz))
        target = hrefs[block_end]

        halls[block_end] = check
        omegas[block_start] = omega
        tplots[block_end] = mag_times[block_end]

        np.save('./hdata/halls%03d.npy' % run, halls)
        np.save('./hdata/omegas%03d.npy' % run, omegas)
        np.save('./hdata/tplots%03d.npy' % run, tplots)

        print('Check, target', check, target)
        #THIS ASSUMES A ZERO INTERCEPT RATHER THAN MINIMISING ANYTHING

        xs.append(omega); ys.append(check - target)

        if abs(check - target) < 1e0:   #Close enough
            go = False
            print('GOOD ENOUGH',  omega, check-target, xs, ys)

        else:   #not close enough, make a better estimate
            print('NOT GOOD ENOUGH', omega, check-target, xs, ys)
            if len(xs) == 1:
                if check > target and omega > minomega:
                    omega = omega - omega_range/4
                elif check < target and omega < maxomega:
                    omega = omega + omega_range/4
                else:
                    go = False
            else:
                #Check for anomalies...
                if (ys[-1] - ys[-2]) / (xs[-1] - xs[-2]) < 0:
                    go = False
                else:
                    nr_target = xs[-1] - ys[-1]*((xs[-1] - xs[-2])/(ys[-1] - ys[-2]))
                    if nr_target > maxomega:
                        omega = np.mean(omegas[block_start-back:block_start])
                        stopnext = True
                        print('OMEGA BIG', nr_target)
                    elif nr_target < minomega:
                        omega = minomega
                        omega = np.mean(omegas[block_start-back:block_start])
                        stopnext = True
                        print('OMEGA SMALL', nr_target)

                    else:
                        print('TARGET',nr_target)
                        omega = nr_target

