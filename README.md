Code for running magnetofrictional simulations using dynamic boundary data, with the new helicity-matching algorithm.

It can be run using either SHARP magnetogram data or synthetic boundary data -- some of which is stored here in the folder 'magnetograms' for use as a test case.

The script 'run.py' is a python wrapper which should contain all you need to effectively set up and run multiple simulations, but the brunt of the calculations are done in fortran using the scripts in ./src/

Steps to install and run:

Clone this repo into your place of choice:

```
git clone 
```
