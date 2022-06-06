---
layout: sidepage
permalink: /summary/
---

# Summary of the library

Here is a summary of the most important functionalities included in TSAOpy. It's advised to read this and then go the notebooks to see examples of the
usage.

## Outline of how TSAOpy works

The first thing needed is the data. The data provided to `tsaopy` would be a set of 1D `numpy` arrays with the $t$, $x(t)$, and (optionally) $v(t)$ measurements. The values must be uniformly distributed in time.

Next, it's necessary to decide what will be the ODE. For each ODE term we will define a `tsaopy` parameter object (details on how to do it later). We also have to define two `tsaopy` parameters for $x_0$ and $v_0$, and optionally we can define another parateter for $X_{eq}$ if the system's equilibrium position is not zero. 

With the data and the parameters it's possible to create the `tsaopy` model object. This object gathers everything needed to set up things like the numerical integrator, the `emcee` methods and functions, and specifically the `emcee` sampler object to run the MCMC chain. 

Once the `tsaopy` model object is built, the method to set up the MCMC sampler is ready and can be run to start the MCMC chain.

After the MCMC chain finishes running then

- If the chain converged it's possible to extract the raw results and process them into the information we were looking for.
- If the chain did not converge we can analyze what happened and run more chains with possibly different set ups and see where we go until we get a converging chain.

From the results of a converged chain one can obtain valuable information such as posterior distributions for each parameter and correlations between each pair of parameters, which is in general the goal of the fitting.

So the overall to do list looks like this

1. Preprocess and import the data.
2. Define `tsaopy` parameters.
3. Define the `tsaopy` model. 
4. Set up and run the MCMC sampler.
5. Extract and process the results.

More details are given in the following sections. 

## TSAOpy parameters

`tsaopy` has two parameter classes, `Fixed` and `Fitting`. The fixed class will represent parameters in the model which the correct value is assumed to be known and therefore considered in the calculations. The fitting class will instead represent parameters for which we don't know the true value. The MCMC method will try to find a probability density

`tsaopy` comes with built in classes to easily set up callable PDFs
