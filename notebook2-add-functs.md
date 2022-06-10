---
layout: sidepage
permalink: /add-functs/
---

Interactive version of the notebook at Google Colab [here](https://colab.research.google.com/drive/1VNKnvGW3U35vLyiSXX9e4blLxyCp4rCD?usp=sharing).

# Notebook 2: Additional functionalities

In this notebook we will show additional functionalities that were not covered in the basic usage notebook. These may be helpful in some cases but generally are not necessary.

We will be using the same model than in the previous notebook.

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

# simulation
a1,b1 = 0.3,1.0
deriv = lambda X : np.array([  X[1],  -a1*X[1]-b1*X[0]  ])

X0 = np.array([1.0,0.0])

result = solve_ivp(deriv,X0,(0,30),0.01)

t,x,v = result[:,0],result[:,1],result[:,2]

# fix time series length
n_data = len(t)
n_out = 1000
step = n_data//n_out

x, v, t = x[::step], v[::step], t[::step]
while len(t) > n_out:
    x, v, t = x[:-1], v[:-1], t[:-1]

# aux noise function
np.random.seed(205)
u_noise, n_noise = 2, 1
noise = lambda : np.random.uniform(-u_noise,u_noise) +\
    np.random.normal(0,n_noise)

# add noise
for i in range(n_out):
    x[i] = x[i] + noise() * 0.05
    v[i] = v[i] + noise() * 0.05
```

## Fitting to $x(t)$ data with a non zero equilibrium point

Suppose our data is shifted such that

$$ x(t) = X(t) - X_{eq} $$

where $x(t)$ is the displacement from the equilibrium position and $X(t)$ is the actual measurement. 

In this case we set up the `tsaopy` model adding a `ptype='ep'` parameter

```
import tsaopy

# import data

t_data,x_data = t,x + 2
x_data_unc = 0.2

# define priors

x0_prior = tsaopy.tools.uniform_prior(0.7,1.3)
v0_prior = tsaopy.tools.normal_prior(0,10)
ep_prior = tsaopy.tools.uniform_prior(1.5,2.5)
a1_prior = tsaopy.tools.normal_prior(0,10)
b1_prior = tsaopy.tools.normal_prior(0,10)

# define tsaopy parameters

x0 = tsaopy.parameters.Fitting(1,'x0',x0_prior)
v0 = tsaopy.parameters.Fitting(0,'v0',v0_prior)
ep = tsaopy.parameters.Fitting(2,'ep',ep_prior)
a1 = tsaopy.parameters.Fitting(0,'a',a1_prior,1)
b1 = tsaopy.parameters.Fitting(0,'b',b1_prior,1)

parameters = [x0,v0,ep,a1,b1]

# define tsaopy model

mymodel = tsaopy.models.PModel(parameters,t_data,x_data,x_data_unc)
sampler = mymodel.setup_sampler(100,300,300)
samples,flat_samples = sampler.get_chain(), sampler.get_chain(flat=True)
labels = mymodel.params_labels
tsaopy.tools.cornerplots(flat_samples,labels)
```

![pic1](https://user-images.githubusercontent.com/94293518/173155321-01fc26a5-d3da-4ca8-b62f-d3541871b51c.png)

## External optimization of the initial fitting values

By applying MCMC we fit the model to the highest likelihood parameter values. However, by using an external optimization we may get a quick estimate that when supplied to the MCMC chain will make it converge much faster. 

Two useful external optimizers are `dual_annealing` and `differential_evolution` from `scipy.optimize`.

In order to do this we first set up the `tsaopy` model as usual.

```
import tsaopy

# import data

t_data,x_data = t,x
x_data_unc = 0.2

# define priors

x0_prior = tsaopy.tools.uniform_prior(0.7,1.3)
v0_prior = tsaopy.tools.normal_prior(0,10)
a1_prior = tsaopy.tools.normal_prior(0,10)
b1_prior = tsaopy.tools.normal_prior(0,10)

# define tsaopy parameters

x0 = tsaopy.parameters.Fitting(1,'x0',x0_prior)
v0 = tsaopy.parameters.Fitting(0,'v0',v0_prior)
a1 = tsaopy.parameters.Fitting(0,'a',a1_prior,1)
b1 = tsaopy.parameters.Fitting(0,'b',b1_prior,1)

parameters = [x0,v0,a1,b1]

# define tsaopy model

mymodel = tsaopy.models.PModel(parameters,t_data,x_data,x_data_unc)
```

Now, we extract the negative log likelihood from the model object using the `neg_ll` method

```
neg_ll = mymodel.neg_ll
```
In this case we will be optimizing `neg_ll` using `dual_annealing`

```
from scipy.optimize import dual_annealing

bounds = [(0.7,1.3),(-10,10),(-10,10),(-10,10)]

da_sol = dual_annealing(neg_ll, bounds)
```

Now we will use the result to set the initial values for an MCMC chain

```
mymodel.update_initvals(da_sol.x)
sampler = mymodel.setup_sampler(100,20,300)
samples,flat_samples = sampler.get_chain(), sampler.get_chain(flat=True)
labels = mymodel.params_labels
tsaopy.tools.traceplots(samples,labels)
```
![pic2](https://user-images.githubusercontent.com/94293518/173155550-b5742802-2502-49f9-b435-eb9b372c5b09.png)

And now with a set of initial values much closer to the posterior means, the chain takes about 50 steps to approach convergence vs 250 in the previous case.
