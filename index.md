# Welcome to TSAOpy

A Python library developed to fit user defined differential equations (with the form of anharmonic oscillators) to time series data that show a roughly periodic behaviour. 

## What it does

Let's assume we have a set of $(t,x(t))$ points making up a time series such as the following

<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/ex_timeseries.png" width="400">

This library will allow you to model the dynamics of the $x(t)$ function as satisfying a differential equation of the form

$$ \ddot{x} + \sum_n a_n \dot{x}|\dot{x}|^{n-1} + \sum_m b_m x^m + \sum_{ij} c_{ij} x^i\dot{x}^j = F_0 \sin{(\omega t + \phi)} $$

along with initial conditions $x(t=0)=x_0$ and $\dot{x}(t=0)=v_0$ required to solve the ODE numerically. Once you define your model by choosing which terms  in the ODE you will consider, the program will fit the model to the data finding the most likely values for each parameter (including initial conditions). This is done using the MCMC method, implemented on the `emcee` library. 

A more detailed explanation can be found at [method page](https://tsaopy.github.io/method/).

## Why doing this?

Running this analysis will allow you to find an ODE that your time series roughly obeys. Some interesting things that you can do with the ODE are
 
1. Modelling of damping forces or potentials. One may associate each term of the ODE with a certain term of the polynomial expansion of a potential or damping force, and thus getting approximations of the behaviour of your system's potential or drag effects.
2. Non linear dynamics analysis. Finding the ODE (or an approximation of it) that your system obeys allows you to use some theoretical tools such as plotting phase portraits and phase space trajectories.


## Set Up guide

In order to make `tsaopy` work you should first be using a Linux PC, specifically it's being developed and tested on Ubuntu 20+ systems. Windows won't work, and we haven't tried it on macOS.

The next step is to make sure `gfortran` is installed in your PC. If it's not, then `sudo apt-get install gfortran` should do the trick. [^1]

You will need some basic Python dependencies like `sys`, `math`, `numpy`, and `matplotlib`. You will also need the more specific dependencies `multiprocessing`, `emcee`, and `corner`. You may install this in the Python console using `pip install` or with an enviroment manager such as Anaconda.

After all this is done you may download the library's files, and then run the 'f2py' script. This will build the module used by the backend. It should take a few seconds and if everything worked a new file should have been created in the working directory. 

At this point `tsaopy` should be working, and you may head to the 'maketestdata' file to create mock time series, or pick up your own data, and then go to [basic usage](https://tsaopy.github.io/basic-usage/) to check how everything works and fit some models. 

### Note: the solution and data files must be in the same directory as the backend files

If your want to keep the notebook files in a lower level directory as provided in the repository, then in the solution file replace

```
import backend as bend
```

for
```
import sys
sys.path.append('..')
import backend as bend
```

## Referencing TSAOpy

We have a Zenodo DOI set up for the repository for the time being

[![DOI](https://zenodo.org/badge/427913804.svg)](https://zenodo.org/badge/latestdoi/427913804)

And this is [my personal ORCID](https://orcid.org/0000-0002-1007-8229)

A possible Bibtex citation is

```
@misc{tsaopy,
      doi = {10.5281/ZENODO.6569849},
      url = {https://tsaopy.github.io/},
      author = {{Scozziero, Sofía A.},
      title = {TSAOpy},
      publisher = {Zenodo},
      year = {2022},
      copyright = {Open Access}}
```

I will probably write an article or something of the sort in the near future. 

[^1]: while the end user should be only concerned with Python, there is some Fortran code used by the backend, and hence why this is necessary.
