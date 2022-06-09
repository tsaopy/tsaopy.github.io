---
layout: sidepage
permalink: /basic-usage/
---

Interactive version of the notebook at Google Colab [here](https://colab.research.google.com/drive/1Ed-LWcUAaPfMktbDX57c9tdV_v8XbyE3?usp=sharing).

# Notebook 1: Basic Usage

In this first notebook we will define a simple model, make a simulation with it, add some noise to the data, and see how to recover the original model using `tsaopy`.

## Making some data

We want to make some data by simulating the following model

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

def solve_ivp(f,x0,t0tf,dt):
    P = x0
    t,tf = t0tf
    x,v = P
    result = [[t,x,v]]
    while t < tf:
        P = rk4(f,P,dt)
        t = t + dt
        x,v = P
        result.append([t,x,v])
    return np.array(result)
```

Then we can run the simulation as

```
# simulation
a1,b1 = 0.3,1.0
deriv = lambda X : np.array([  X[1],  -a1*X[1]-b1*X[0]  ])

X0 = np.array([1.0,0.0])

result = solve_ivp(deriv,X0,(0,30),0.01)

t,x,v = result[:,0],result[:,1],result[:,2]
```
We'll add another block to fix the length of the arrays to a given number we want

```
# fix time series length
n_data = len(t)
n_out = 1000
step = n_data//n_out

x, v, t = x[::step], v[::step], t[::step]
while len(t) > n_out:
    x, v, t = x[:-1], v[:-1], t[:-1]
```

Now, let's plot what we have so far

```
import matplotlib.pyplot as plt

# plot
plt.figure(figsize=(7, 5), dpi=100)
plt.plot(t, x, color = 'tab:red', label='position')
plt.plot(t, v, color='tab:purple', label='velocity')
plt.legend()
plt.show()
```
![pic1](https://user-images.githubusercontent.com/94293518/172937595-3294ca14-e1d2-4b06-8dac-9d781dfebc9c.png)

Next we will add noise to the data and plot it

```
# aux noise function
np.random.seed(205)
u_noise, n_noise = 2, 1
noise = lambda : np.random.uniform(-u_noise,u_noise) +\
    np.random.normal(0,n_noise)

# add noise
for i in range(n_out):
    x[i] = x[i] + noise() * 0.05
    v[i] = v[i] + noise() * 0.05
    
# plot
plt.figure(figsize=(7, 5), dpi=100)
plt.plot(t, x, color = 'tab:red', label='position')
plt.plot(t, v, color='tab:purple', label='velocity')
plt.legend()
plt.show()
```
![pic2](https://user-images.githubusercontent.com/94293518/172937829-cc86c714-5fb7-4336-858c-018689d7d6ac.png)

## Fitting the model with `tsaopy`

Now that we have our data, we will try to recover the original equation.

The first thing needed to set up our `tsaopy` model is importing the data, including both the measurements and the uncertainty.

```
# import data

# if we haven't generated the data here
# we would be importing it from the source

t_data,x_data = t,x
```
We need to define the uncertainty of the $x(t)$ measurements. It can be either a single number representing the uncertainty of all measurements or an array of the same length as `x_data` with an uncertainty value for each value of $x$. In this case we will use a baseline for all masurements.

```
# define uncertainty

x_data_unc = 0.2
```

Next thing we need is to define the parameters for the model. Defining `tsaopy` parameter objects will allow `tsaopy` to build the ODE and everything needed for the fitting. The parameters we are considering for this model are $x_0$, $v_0$ (since we are solving an ODE initial conditions are always required), and the coefficients from the ODE $a_1$ and $b_1$. 

MCMC needs priors for each parameter, we can set them up easily with the prior classes defined in `tsaopy.tools`. For $x_0$ we will assume its value is between $0.7$ and $1.3$ and use a uniform prior. For $v_0$, $a_1$, and $b_1$ we will define gaussian priors centered at 0.

```
# define priors

x0_prior = tsaopy.tools.uniform_prior(0.7,1.3)
v0_prior = tsaopy.tools.normal_prior(0,10)
a1_prior = tsaopy.tools.normal_prior(0,10)
b1_prior = tsaopy.tools.normal_prior(0,10)
```

With the priors defined we can define `tsaopy` parameters. 

```
# define tsaopy parameters

x0 = tsaopy.parameters.Fitting(1,'x0',x0_prior)
v0 = tsaopy.parameters.Fitting(0,'v0',v0_prior)
a1 = tsaopy.parameters.Fitting(0,'a',a1_prior,1)
b1 = tsaopy.parameters.Fitting(0,'b',b1_prior,1)

parameters = [x0,v0,a1,b1]
```

With the parameters and the data we can set up the `tsaopy` model object as

```
mymodel = tsaopy.models.PModel(parameters,t_data,x_data,x_data_unc)
```

Now we can start running MCMC chains with it. Let's start by running a test chain, we will use 100 walkers, and will have a shorter burn in, so we can see how the sampler explores the space. 

```
sampler = mymodel.setup_sampler(100,20,300)
```

Now let's extract the results and make traceplots

```
samples,flat_samples = sampler.get_chain(), sampler.get_chain(flat=True)
labels = mymodel.params_labels
tsaopy.tools.traceplots(samples,labels)
```
![pic3](https://user-images.githubusercontent.com/94293518/172950111-79208438-852f-4976-9565-b0704967b637.png)

We can see from the traceplot that the chain looks like it has converged in about 300 steps. Let's do another run with a 300 steps burn in, so we can store samples for the result in the production phase.

```
sampler = mymodel.setup_sampler(100,300,300)
samples,flat_samples = sampler.get_chain(), sampler.get_chain(flat=True)
tsaopy.tools.traceplots(samples,labels)
```

![pic4](https://user-images.githubusercontent.com/94293518/172951087-4a79c4e6-1ace-4f47-9928-f4e6397ef79f.png)

New traceplots suggest that the samples may belong to a converged chain. Let's plot the autocorrelation functions to double check

```
tsaopy.tools.autocplots(flat_samples,labels)
```
![pic5](https://user-images.githubusercontent.com/94293518/172951503-ce32945a-a1ac-4fa2-9eb9-fa558477c929.png)


On the other hand, autocorrelation plots suggest that the chain may be too short, so we run another chain with a longer production phase

```
sampler = mymodel.setup_sampler(100,300,1000)
samples,flat_samples = sampler.get_chain(), sampler.get_chain(flat=True)
tsaopy.tools.traceplots(samples,labels)
tsaopy.tools.autocplots(flat_samples,labels)
```
![pic6](https://user-images.githubusercontent.com/94293518/172951684-afb8b098-27c9-4901-9821-dfbb7b77d005.png)
![pic7](https://user-images.githubusercontent.com/94293518/172951708-5982b02c-0144-4a38-95d9-859ab24fe0a6.png)


And after convincing ourselves that the chain has converged we show the results in cornerplots

![pic8](https://user-images.githubusercontent.com/94293518/172951737-ab98334b-1717-4b4a-b4ea-7a62f874e271.png)
