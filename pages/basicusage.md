---
layout: sidepage
permalink: /basic-usage/
---

Interactive version of the notebook at Google Colab [here](https://colab.research.google.com/drive/1Ed-LWcUAaPfMktbDX57c9tdV_v8XbyE3?usp=sharing).

# Basic usage

In this notebook we will define a simple model, make a simulation with it, add some noise to the data, and see how to recover the original model using `tsaopy`.

# Making some data

We want to make some data by solving the ODE

$$ \ddot{x} + a_1\dot{x} + b_1x = 0\quad; \qquad x_0=1 \qquad v_0=0 \qquad a_1=0.3 \qquad b_1=1 $$

We will start by setting up the numerical integrator using Runge-Kutta 4. A simple implementation is

```
import numpy as np

# numerical solver
def rk4(f,x,dx):
    k1 = dx*f(x)
    k2 = dx*f(x+0.5*k1)
    k3 = dx*f(x+0.5*k2)
    k4 = dx*f(x+k3)
    return x + (k1+k4+2*k2+2*k3)/6

def solve_ivp(f, x0, t0tf, dt):
    t, tf = t0tf
    x, v = x0
    result = [[t, x, v]]
    while t < tf:
        x, v = rk4(f, (x, v), dt)
        t = t + dt
        result.append([t, x, v])
    return np.array(result)
```

Then we run the simulation as

```
# simulation
a1, b1 = .3, 1.0
deriv = lambda X : np.array([X[1],
                            -a1*X[1]-b1*X[0]])

X0 = np.array([1.0,0.0])

result = solve_ivp(deriv, X0, (0, 30.0), .01)

t,x,v = result[:,0],result[:,1],result[:,2]
```
We add another block to fix the length of the arrays to a given number we want

```
# fix time series length
n_data = len(t)
n_out = 1000
step = n_data // n_out

x, v, t = x[::step], v[::step], t[::step]
while len(t) > n_out:
    x, v, t = x[:-1], v[:-1], t[:-1]
```

We plot what we have so far

```
import matplotlib.pyplot as plt

# plot
plt.figure(dpi=100)
plt.plot(t, x, color = 'tab:red', label='position')
plt.plot(t, v, color='tab:purple', label='velocity')
plt.legend()
plt.show()
```
![index](https://user-images.githubusercontent.com/94293518/180663961-3b2d5902-8a72-4855-b30f-aa7209ae205c.png)

Next we add noise to the data and plot it

```
# aux noise function
np.random.seed(205)
u_noise, n_noise = 2, 1
noise = lambda: np.random.uniform(-u_noise, u_noise) +\
    np.random.normal(0, n_noise)

# add noise
for i in range(n_out):
    x[i] = x[i] + noise() * .05
    v[i] = v[i] + noise() * .05
    
# plot
plt.figure(dpi=100)
plt.plot(t, x, color = 'tab:red', label='position')
plt.plot(t, v, color='tab:purple', label='velocity')
plt.legend()
plt.show()
```
![index](https://user-images.githubusercontent.com/94293518/180663985-6e819aa4-8158-4b88-a2ae-a7e2d52e7084.png)

# Fitting the model with TSAOpy

Now that we have our mock data, we will try to recover the original equation. We are assuming the model we used earlier to generate the data is correct, but now we pretend to not know the values of the parameters $x_0$, $v_0$, $a_1$, and $b_1$, an we want to find them by fitting the model to the data.

## Importing data

The first thing needed to set up our `tsaopy` model is importing the data, including both the measurements and their uncertainty.

```
# import data

# if we haven't generated the data here
# we import it from the source

t_data, x_data = t, x
```
To define the uncertainty of the measurements, it can be either a single number representing the uncertainty of all measurements or an array of the same shape as `x_data` with an uncertainty value for each value of x in the array. In this case we will use a single number for all masurements.

```
# define uncertainty

x_data_unc = .2
```

## Set up some priors

The next thing we need is to set up priors for the parameters to fit.

We will start by setting up the prior for $x_0$ and $v_0$. Since we are pretending to not know anything about their values, we will just use normal distributions centered at $0$ with a large SD.

```
# define priors

import quickemcee as qmc

x0_prior = qmc.utils.normal_prior(.0, 10.0)
v0_prior = qmc.utils.normal_prior(.0, 10.0)
```

## Set up TSAOpy objects

With the data and the priors we set up a `tsaopy` 'Event' object, which summarizes the information about the measurements and the initial conditions.

```
# define tsaopy event

import tsaopy

# this dictionary takes the variable names x0 and v0 as keys
# and each key has the prior for that parameter as its value
event_dict = {'x0': x0_prior,
              'v0': v0_prior}

event = tsaopy.events.Event(params=event_dict, t_data=t_data, x_data=x_data,
                            x_sigma=x_data_unc)
```

Similarly, we now make a dictionary for the ODE coefficients, but now the values for each key will be a list which contains a touple for each coefficient. The touple for each coefficient has its index as first element and its prior as second element.

```
a1_prior = qmc.utils.normal_prior(.0, 10.0)
b1_prior = qmc.utils.normal_prior(.0, 10.0)

ode_dict = {'a': [(1, a1_prior)],
            'b': [(1, b1_prior)]}
```

Now we can build the `tsaopy` 'Model' object, which will automatically set up everything for running MCMC. It takes as arguments the `ode_dict` dictionary, and a list containing `tsaopy` 'Event' objects. The purpose is to allow the user to fit different sets of measurements that may have different initial conditions, equilibrium points, etc. 

In this case we will only use the one we defined earlier however.

```
tsaopymodel = tsaopy.models.Model(ode_coefs=ode_dict, events=[event])
```

## Optimize the initial values

Before running MCMC, we will optimize the initial values of the chain by optimizing the negative logarithmic likelihood function.

```
def neg_ll(coords):
  return - tsaopymodel._log_likelihood(coords)

from scipy.optimize import minimize

sol = minimize(neg_ll, x0=np.zeros(4))

if sol.success:
  x0 = sol.x
else:
  print("The external optimizer didn't converge.")
```

## Run MCMC chains

Now we can start running MCMC chains by setting up a `quickemcee` 'Model' object with it(this is done in one line, don't worry). Let's start running a test chain, we will use 100 walkers, and a short burn in so we can see how the sampler explores the space. 

```
mcmcmodel = tsaopymodel.setup_mcmc_model()
mcmcsampler = mcmcmodel.run_chain(100, 10, 100,
                                  init_x=x0)
```
```
>>>Running burn-in...
>>>100%|██████████| 10/10 [00:00<00:00, 24.77it/s]
>>>Running production...
>>>100%|██████████| 100/100 [00:07<00:00, 13.04it/s]
```

## Analyze results

Now let's extract the results

```
samples,flat_samples = mcmcsampler.get_chain(), mcmcsampler.get_chain(flat=True)
labels = tsaopymodel.paramslabels
```

And make traceplots

```
# traceplots
qmc.utils.traceplots(samples, labels)
```
![index](https://user-images.githubusercontent.com/94293518/180664284-61a01bf0-6095-442d-aa4f-2312cce58495.png)

We see that a burn in of about 100 steps might be enough for the chain to converge, so we run another chain.

```
mcmcmodel = tsaopymodel.setup_mcmc_model()
mcmcsampler = mcmcmodel.run_chain(100, 100, 100,
                                  init_x=x0)
samples,flat_samples = mcmcsampler.get_chain(), mcmcsampler.get_chain(flat=True)
```
```
>>>Running burn-in...
>>>100%|██████████| 100/100 [00:10<00:00,  9.41it/s]
>>>Running production...
>>>100%|██████████| 100/100 [00:03<00:00, 25.31it/s]
```
And re do the traceplots
```
# traceplots
qmc.utils.traceplots(samples, labels)
```
![index](https://user-images.githubusercontent.com/94293518/180664333-f932b9e4-19fe-419c-8c58-82518533cc10.png)

The chain looks converged now, so we make autocorrelation plots.

```
# autocplots
qmc.utils.autocplots(flat_samples, labels)
```
![index](https://user-images.githubusercontent.com/94293518/180664349-2fba4d42-0c44-4e0e-a00e-f7b7aff9f45e.png)

Here we see that we could use more samples, so we re do the chain with a longer production phase and re do the autocorrelation plots.

```
mcmcmodel = tsaopymodel.setup_mcmc_model()
mcmcsampler = mcmcmodel.run_chain(100, 100, 500,
                                  init_x=x0, workers=2)
samples,flat_samples = mcmcsampler.get_chain(), mcmcsampler.get_chain(flat=True)

# autocplots
qmc.utils.autocplots(flat_samples, labels)
```
```
>>>Running burn-in...
>>>100%|██████████| 100/100 [00:04<00:00, 21.78it/s]
>>>Running production...
>>>100%|██████████| 500/500 [00:20<00:00, 23.87it/s]
```
![index](https://user-images.githubusercontent.com/94293518/180664382-1aea0e43-5414-4696-a056-447970e537ed.png)

Finally we show the resulting posteriors in cornerplots.

```
# cornerplots
qmc.utils.cornerplots(flat_samples, labels)
```
![index](https://user-images.githubusercontent.com/94293518/180664395-6e83f29a-903a-4eac-963f-ec33d2048116.png)
