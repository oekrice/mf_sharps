Code for running magnetofrictional simulations using dynamic boundary data, with the new helicity-matching algorithm.

It can be run using either SHARP magnetogram data or synthetic boundary data -- some of which is stored here in the folder 'magnetograms' for use as a test case.

The script 'run.py' is a python wrapper which should contain all you need to effectively set up and run multiple simulations, but the brunt of the calculations are done in fortran using the scripts in ./src/

## Steps to install and run:

# 1. Clone this repo into your place of choice:

```
git clone https://github.com/oekrice/mf_sharps.git
```

which should download everything. I'd recommend doing this rather than copying directly as there are lots of files dotted around. This will also download the sample magnetograms (around 100MB I think). 

# 2. Check the Makefile works by running:

```
make clean
make
```

If not, then you'll need to open the Makefile and change the lines at the top for wherever your respective MPI and netcdf files are located (currently set up for local maths machine).

# 3. Change the location of the temporary and output files

These are currently set up to save to /extra/tmp/ . Lots of space is required here. 

You'll need to change lines 54 to 56 (ish) in 'run.py', along with the location of your mpi executable, and the 'data_directory_root' on line 31 in src/main.f90. I *think* that's all of them...  

All other temporary files will be saved locally in the repo folder, but these aren't very big.

# 4. Set up simulation 

The (rather large) preamble of 'run.py' should contain everything, and it's commented such that hopefully it makes sense. Choices to use the synthetic or SHARP data, region or sharp ids etc. are here. Of note are 'adapt_omega' which turn on/off the helicity matching steps, and 'dothings' which recalculates all the initial and boundary stuff. If you're doing multiple runs with the same boundary setup, set this to false. 

# 5. Input system parameters

These are from line 74 onwards and hopefully are familiar. Shouldn't need to touch below line 160...

# 6. Run things

Run using
```
python run.py 0 1
```
to test. The first argument is a run id, allowing for multiple simulations saving to diferent places (these are called in the diagnostics etc. and plotting tools) and the second argument is the number of cores to run the fortran with. Sensible even numbers should all work.

Outputs will be saved as netcdfs to the file specified above. Basic diagnostics can be done with the scripts

```
python helicity_adapt.py
python diagnostics_new.py
```

though within these scripts you'll need to change the run ids if it's not zero.

# 7. Field line tracing

This uses the other repo fltrace, which I've kept separate from this one for good reasons. The files are completely compatible though.

Enjoy.

