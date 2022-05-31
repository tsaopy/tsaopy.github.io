---
layout: sidepage
permalink: /method/
---

# Method

## The differential equation

Picking up from where we left, we assumed that we have a set of $(t,x(t))$ points that make up a time series, like this one which we showed earlier

<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/ex_timeseries.png" width="700">

We also proposed that the $x(t)$ function satisfies a differential equation of the form

$$ \ddot{x} + \sum_n a_n \dot{x}|\dot{x}|^{n-1} + \sum_m b_m x^m + \sum_{ij} c_{ij} x^i\dot{x}^j = F_0 \sin{(\omega t + \phi)} $$

subject to the initial conditions $x(t=0)=x_0$ and $v(t=0)=v_0$. One think to also keep in mind is that in the ODE $x(t)$ is not the actual position of the oscillator, but rather it's the distance to the equilibrium position. So if the equilibrium position is not equal to zero then this approach will fail. In order to fix that we assume that we can switch between **position** from the lab frame of reference and **distance to the equilibrium position**, as if you set the zero at the equilibrium point, as $X(t) = X_{eq} + x(t)$. For all intents and purposes assume that the function that goes into the ODE is $x(t)$ and we get $X(t)$ back by adding $X_{eq}$.

Now we'll make some more remarks about the models.

The terms proportional to $x$ can be associated with the potential energy, we may think of them as a polynomial expansion of the force caused by the potential, as if the potential has a series expansion and we truncated it at some degree.

The same goes for the $\dot{x}$ terms which we associate to a damping or drag force. Notice that these terms are written as $\dot{x}\mid\dot{x}\mid ^{n-1}$ instead of $\dot{x}^n$. The purpose of this tweak is that the absolute value of the damping force will always be proportional to some power of the absolute value of the velocity, and will have a direction opposite to the velocity, as we'd tipically expect from a drag force. 

We included the possibility of adding a sinusoidal driving force to the model. It is not possible to add driving forces with a different law at the moment. 

It's also not possible to set up relations between each parameter, such as $b_1 = 4 a_1$. Each parameter will be treated as independent from the others (you may however, find correlations after getting the results of the fitting). 


## The MCMC fitting

MCMC is a method to find a posterior distribution for a set of parameters, given some previous knowledge about the parameters (in the form of prior probability distributions), and new observations which are subject to the model. An in depth explanation is given in the official [emcee documentation](https://emcee.readthedocs.io/en/stable/tutorials/line/).

So, in order to apply the MCMC method we need three things, a model, a prior probability distribution for each parameter of the model, and a set of observations to which we'll fit the model. 

### More about priors

When approaching MCMC for the first time the concepts of priors and posteriors might not be very clear. When we talk about those two, the mathematical objects which describe them are probability density functions, such as the normal (aka Gaussian) that might be familiar. So both priors and posteriors are PDFs, but the prior is the PDF that represented what we knew about the parameter before applying the MCMC method, and the posterior is the PDF that represents our knowledge about the parameter after MCMC, it's "the result" of solving the MCMC problem.

Here are two examples of setting up priors. Suppose we have an object, and we want to express what we know about its mass. Before weighing it we know the mass of a regular object must be a positive finite number, so:

1. It must be greater than zero.
2. It must be lower than some upper bound $M$.

With that knowledge we can set up a prior PDF to express what we now, specifically it will be a uniform distribution with the form:

$$ p(m) = \begin{cases} c \qquad \text{if }0<m<M \\ 0 \qquad \text{otherwise} \end{cases} $$

### Some specifics about the MCMC sampler

When we set up the MCMC sampler we will have to specify three values, which we call walkers, burn in steps, and production steps.

These are explained in `emcee` docs but, summing up, walkers are the number of chains that we are running at the same time. So, if we have 10 walkers, at each step we are proposing 10 new samples, if we have 50 walkers we propose 50 new samples at each step and so on. 

Now, burn in and production steps are the number of steps that the sampler will do at burn in stage and production stage, respectively. Ideally, during burn in stage the sampler will be drawing samples while exploring the possible values for your parameters until it finds a region where the likelihood is at a maximum, and we are accepting a lot of the newly drawn samples. In the production stage, assuming the chain has converged, we will simply be drawing samples that will belong to the posterior distribution we were looking for, and we will run it for as long as we need to get the number of samples we want. 
