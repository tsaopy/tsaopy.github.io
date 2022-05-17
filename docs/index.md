# Welcome to TSAOpy

A Python library developed to fit user defined differential equations, with the form of anharmonic oscillators, to time series data that show a roughly periodic behaviour.

## Set Up guide

In order to make `tsaopy` work you should first be using a Linux PC, specifically it's being developed and tested on Ubuntu 20+ systems. Windows won't work, and we haven't tried it on macOS.

The next step is to make sure `gfortran` is installed in your PC. If it's not, then `sudo apt-get install gfortran` should do the trick. [^1]

You will need some basic Python dependencies like `sys`, `math`, `numpy`, and `matplotlib`. You will also need the more specific dependencies `multiprocessing`, `emcee`, and `corner`. You may install this in the Python console using `pip install` or with an enviroment manager such as Anaconda.

After all this is done you may download the library's files, and then run the 'f2py' script. This will build the module used by the backend. It should take a few seconds and if everything worked a new file should have been created in the working directory. 

At this point `tsaopy` should be working, and you may head to the 'maketestdata' file to create mock time series, or pick up your own data, and then go to the 'notebook' file to check the basic usage and fit some models. 

[^1]: while the end user should be only concerned with Python, there is some Fortran code used by the backend, and hence why this is necessary.
