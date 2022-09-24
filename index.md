---
layout: homepage
---

# Welcome to TSAOpy's site

Time Series by Anharmonic Oscillators is a Python library developed to fit user defined differential equations (with the form of anharmonic oscillators) to time series data(aimed at series with oscillating behaviour). 

## What it does

Let's assume we have sets of points making up time series such as the following

![index](https://user-images.githubusercontent.com/94293518/180660624-b9c43ffd-2f0f-48b3-96b8-bac9a13a1f91.png)

This library will allow the user to model the dynamics of the $x(t)$ function with a differential equation of the form

$$ \ddot{x} + \sum_n a_n \dot{x}|\dot{x}|^{n-1} + \sum_m b_m x^m + \sum_{ij} c_{ij} x^i\dot{x}^j = F(t) $$

along with initial conditions $x(t=0)=x_0$ and $\dot{x}(t=0)=v_0$ required to solve the ODE numerically. The program will fit the model to the datasets finding the most likely values for each parameters. The fitting is done using a maximum likelihood estimation, and the MCMC method, implemented in the `emcee` library, is used to find a posterior distribution for the parameters. 

## Why doing this?

Running this analysis allows one to find an ODE that a(or a set of) time series roughly obeys. Some interesting things that can be done with the ODE are
 
1. Modelling of damping forces or potentials. One may associate each term of the ODE with a certain term of the polynomial expansion of a potential or damping force, and thus getting approximations of the behaviour of your system's potential or drag effects.
2. Non linear dynamics analysis. Finding the ODE that the system obeys allows using some theoretical tools such as plotting phase portraits and phase space trajectories, finding limit  cylces, analyzing stability, energy conservation, etc.


## Set Up

- In order to make `tsaopy` work one should first be using a Linux PC, specifically it's being developed and tested on Ubuntu 20+ systems. Windows won't work, and we haven't tried it on macOS.

- Necessary Python dependencies are `numpy`, and another package I developed(largely to serve as a base for `tsaopy`) which is `quickemcee`, which also depends on `scipy`, `matplotlib`, `emcee`, and `corner`. 

- `tsaopy` uses a Fortran submodule in its backend. In order to build it, it's necessary that there's a Fortran compiler, such as `gfortran`, installed and properly pathed. If there isn't, then `sudo apt-get install gfortran` should install a compiler. 

- With this now it's possible to install `tsaopy` using `pip install tsaopy`. 

At this point `tsaopy` should be installed.

### Note: some users may run into trouble if path variables for certain numpy submodules are not properly set up, and I can't help with that. However, it should work out of the box if you install things in a new conda enviroment, and all the dependencies including Python itself are up to date. 

If you have any problems during installation please make an issue in the Github repo with all the information you can gather.

### Test installation

After running `pip install tsaopy` try opening a Python console and run `import tsaopy`. It may take a few seconds (the backend submodules are compiled the first time you import `tsaopy`). Should any errors arise, please report the issue.

If you can import `tsaopy` succesfully try running the basic test script. Go to the project's repository test folder and either

- Download 'basictest.py' and run it from the Linux terminal with `python basictest.py`, or
- Copy the contents of the file and run them on a Python console.

If everything is working properly something like this should be displayed

![index](https://user-images.githubusercontent.com/94293518/192122736-f272aa9b-5106-41be-bb1c-b4c560181531.png)


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

There will probably be an article or something of the sorts in the near future. 
