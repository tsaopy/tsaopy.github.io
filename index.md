# Welcome to TSAOpy

A Python library developed to fit user defined differential equations, with the form of anharmonic oscillators, to time series data that show a roughly periodic behaviour. 

## What does it do?

Let's assume we have a set of $(t,x(t))$ points making up a time series such as the following

<img src="https://github.com/tsaopy/tsaopy.github.io/blob/main/assets/ex_timeseries.png" width="200">

This library will allow you to model the dynamics of the $x(t)$ function as satisfying a differential equation of the form

$$ \ddot{x} + \sum_n a_n \dot{x}|\dot{x}|^{n-1} + \sum_m b_m x^m + \sum_{ij} c_{ij} x^i\dot{x}^j = F_0 \sin{(\omega t + \phi)} $$

along with initial conditions $x(t=0)=x_0$ and $\dot{x}(t=0)=v_0$ required to solve the ODE numerically. Once you define your model by choosing which terms you will consider, the program will fit the model to the data finding the most likely values for each parameter (including initial conditions). This is done using the MCMC method, implemented on the `emcee` library. 

A more detailed explanation can be found at method page.

## Set Up guide

In order to make `tsaopy` work you should first be using a Linux PC, specifically it's being developed and tested on Ubuntu 20+ systems. Windows won't work, and we haven't tried it on macOS.

The next step is to make sure `gfortran` is installed in your PC. If it's not, then `sudo apt-get install gfortran` should do the trick. [^1]

You will need some basic Python dependencies like `sys`, `math`, `numpy`, and `matplotlib`. You will also need the more specific dependencies `multiprocessing`, `emcee`, and `corner`. You may install this in the Python console using `pip install` or with an enviroment manager such as Anaconda.

After all this is done you may download the library's files, and then run the 'f2py' script. This will build the module used by the backend. It should take a few seconds and if everything worked a new file should have been created in the working directory. 

At this point `tsaopy` should be working, and you may head to the 'maketestdata' file to create mock time series, or pick up your own data, and then go to the 'notebook' file to check the basic usage and fit some models. 

[^1]: while the end user should be only concerned with Python, there is some Fortran code used by the backend, and hence why this is necessary.
