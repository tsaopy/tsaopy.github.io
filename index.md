---
layout: homepage
---

# Welcome to TSAOpy's site

Time Series by Anharmonic Oscillators is a Python library developed to fit user defined differential equations (with the form of anharmonic oscillators) to time series data(aimed at series with oscillating behaviour). 

## What it does

Let's assume we have a set of $(t,x(t))$ points making up a time series such as the following

<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/ex_timeseries.png" width="500">

This library will allow the user to model the dynamics of the $x(t)$ function as satisfying a differential equation of the form

$$ \ddot{x} + \sum_n a_n \dot{x}|\dot{x}|^{n-1} + \sum_m b_m x^m + \sum_{ij} c_{ij} x^i\dot{x}^j = F_0 \sin{(\omega t + \phi)} $$

along with initial conditions $x(t=0)=x_0$ and $\dot{x}(t=0)=v_0$ required to solve the ODE numerically. Once the model is defined by choosing which terms  in the ODE will be considered, the program will fit the model to the datasets finding the most likely values for each parameter (including initial conditions). The fitting is done using the MCMC method, implemented in the `emcee` library. 

More details are given in [methodology](https://tsaopy.github.io/methodology/).

## Why doing this?

Running this analysis allows one to find an ODE that a(or a set of) time series roughly obeys. Some interesting things that can be done with the ODE are
 
1. Modelling of damping forces or potentials. One may associate each term of the ODE with a certain term of the polynomial expansion of a potential or damping force, and thus getting approximations of the behaviour of your system's potential or drag effects.
2. Non linear dynamics analysis. Finding the ODE that the system obeys allows using some theoretical tools such as plotting phase portraits and phase space trajectories, finding limit  cylces, analyzing stability, energy conservation, etc.


## Set Up

- In order to make `tsaopy` work one should first be using a Linux PC, specifically it's being developed and tested on Ubuntu 20+ systems. Windows won't work, and we haven't tried it on macOS.

- Necessary Python dependencies are `numpy`, and another package I developed(largely to serve as a base for newer versions of `tsaopy`) which is `quickemcee`, which also depends on `scipy`, `matplotlib`, `emcee`, and `corner`. 

- `tsaopy` uses a Fortran module in its backend. In order to build this module it is necessary that there is a Fortran compiler, such as `gfortran`, installed and properly pathed. If it's not, then `sudo apt-get install gfortran` should install the compiler. 

- With this now it's possible to install `tsaopy` using `pip install tsaopy`. 

At this point `tsaopy` should be installed, and you may head to the [basic usage](https://tsaopy.github.io/basic-usage/) notebook to check how everything works. 

### Note: some users may run into trouble if path variables for certain numpy submodules are not properly defined, and I can't help with that. However, it should work out of the box if you install things in a new conda enviroment, and all the dependencies including Python itself are up to date. 

If you have any problems during installation please make an issue in the Github repo with all the information you can gather.

### Test installation

**Remember to check that you are running the console in the same enviroment that you installed TSAOpy.**

After running `pip install tsaopy` try opening a Python console and run `import tsaopy`. It may take a few seconds (some backend modules are built the first time you import `tsaopy`). Should any errors arise, please report the issue.

If you can import `tsaopy` succesfully try running the basic test script. Go to the project's repository test folder and either

- Download 'basictest.py' and run it from the Linux terminal with `python basictest.py`, or
- Copy the contents of the file and run them on a Python console.

If everything is working properly something like this should be displayed

<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/mainpage_sample_pic.png" width="400">

## Referencing TSAOpy

We have a Zenodo DOI set up for the repository for the time being

[![DOI](https://zenodo.org/badge/427913804.svg)](https://zenodo.org/badge/latestdoi/427913804)

And this is [my personal ORCID](https://orcid.org/0000-0002-1007-8229)

The following Bibtex citation is recommended

```
@software{tsaopy,
  author       = {Scozziero, Sofía Anna},
  title        = {The TSAOpy library},
  month        = may,
  year         = 2022,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.6569848},
  url          = {https://tsaopy.github.io/}
}
```

There will probably be an article or something of the sort in the near future. 
