---
layout: sidepage
permalink: /methodology/
---

# Methodology

## The model

In `tsaopy` we assume that the time series can be modelled as a function $x(t)$ that satisfies a differential equation of the form

$$ \ddot{x} + \sum_n a_n \dot{x}|\dot{x}|^{n-1} + \sum_m b_m x^m + \sum_{ij} c_{ij} x^i\dot{x}^j = F(t, \mathfrak{f}_1, \dotsc, \mathfrak{f}_{N_f}) $$

and we want to find the $a_n$, $b_n$, $c_{ij}$, $\mathfrak{f}_n$ values that best fit our time series data.

## Computing trajectories

By proposing a set of parameter values and initial conditions we can compute a trajectory by numerically solving the differential equation.

* We use Runge-Kutta implemented in Fortran and wrapped to be called from Python with `numpy.f2py`.

## The error function

Once we compute a trajectory with the proposed parameter values ($\theta$), we can compare the simulated data points($f_\theta$) with the real measurements($Y$) with an error or cost function. By default `tsaopy` uses the 'negative logarithmic likelihood' of the measurements given the proposed parameters

$$ - \log{p(Y|\theta)} = \frac{1}{2}\sum_{i}\frac{(f_{\theta, i}-Y_i)^2}{\sigma_i^2} $$

minimizing this function is equivalent to performing a maximum likelihood estimation (MLE). The function can be minimized with any multivariate minimizer such as `minimize` from  `scipy.optimize`.

## Finding the posterior distribution with MCMC

After MLE, we use the MCMC implementation `emcee` to find the posterior distribution of the parameters. This allows us to

* find a value as well as error bars for each parameter
* analyze correlations between parameters
