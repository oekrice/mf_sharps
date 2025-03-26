Code for running magnetofrictional simulations using dynamic boundary data, with the new helicity-matching algorithm.

It can be run using either SHARP magnetogram data or synthetic boundary data -- some of which is stored here in the folder 'magnetograms' for use as a test case.

The script 'run.py' is a python wrapper which should contain all you need to effectively set up and run multiple simulations, but the brunt of the calculations are done in fortran using the scripts in ./src/

Steps to install and run:

Clone this repo into your place of choice:

```
git clone https://github.com/oekrice/mf_sharps.git
```

which should download everything. I'd recommend doing this rather than copying directly as there are lots of files dotted around. This will also download the sample magnetograms (around 100MB I think). 

Check the Makefile works by running:

```
make clean
make
```

If not, then you'll need to open the Makefile and change the lines at the top for wherever your respective MPI and netcdf files are located (currently set up for local maths machine).

