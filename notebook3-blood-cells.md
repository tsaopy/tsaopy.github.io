---
layout: custompage
permalink: /blood-cells/
---

# Notebook 3: Blood cells experiment

Remember that the experiment data file, and other files summing up the work on this notebook can be found at https://github.com/tsaopy/tsaopy/tree/main/notebook3.

## Sketch of the problem

I've been given this and some more datasets by professor Castellini at my home faculty. I believe the experiments were run by colleagues of him, [here's an article reporting some of their work](https://afan.df.uba.ar/journal/index.php/analesafa/article/view/129). From what I've gathered, these are measurements of transverse and longitudinal deformations of red blood cells which were excited by some driving force. Here's the data series that I chose to analyze

<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb3_pic1.png" width="600">

I don't have any more information about this particular set of data points rather than the fact that there was a driving force. So I'm going to begin the analysis with the simplest physical oscillator with a driving force, the ODE will be

$$ \ddot{x} + a_1\dot{x} + b_1x = F_0\sin{(\omega t+\phi)} $$

I will start by assuming that $a_1$ and $b_1$ could be anything, so I will set up normal priors centered at zero with a large SD for those two. Judging from the plot, the oscilation was not near an equilibrium point at the beggining of the time series, so there should be a non zero initial velocity. Also, the oscillator seems to be moving upwards so I'm assuming this velocity is positive. With respect to $x_0$ I'll assume it's negative or slightly positive. So the priors for those parameters will be

```
x0_prior = bend.uniform_prior(-500,50)
v0_prior = bend.uniform_prior(-200.0,2000.0)
a1_prior = bend.normal_prior(0.0,1000.0)
b1_prior = bend.normal_prior(0.0,1000.0)
```

With respect to the driving force, $F_0$ and $\omega$ are strictly positive, and $\phi\in[-\pi,\pi]$, so I'll use uniform priors for those. If you look at the graph the oscillation seems to have period of around $T\approx0.6$, so we can try assuming a frequency around $\omega=2\pi f=2\pi/T\approx10$. So, for $\omega$ we will start the chain at 10, and will be using a uniform prior between 0 and, let's say, 100 so we give it some room to walk. The oscillation seems to have an amplitude of around 500, so $F_0$ could be expected to sit around $5000$ since $A\propto F_0/\omega$, and again let's allow it to walk over an order of magnitude. So, we will start with the next priors for the driving force

```
f_prior = bend.uniform_prior(0.0,50000.0)
w_prior = bend.uniform_prior(0.0,100.0)
p_prior = bend.uniform_prior(-np.pi,np.pi)
```

And then we set up the parameters as

```
x0 = bend.FittingParameter(-100.0,'x0',1,x0_prior)
v0 = bend.FittingParameter(100.0,'v0',1,v0_prior)
a1 = bend.FittingParameter(0.0, 'a', 1, a1_prior)
b1 = bend.FittingParameter(0.0,'b',1,b1_prior)
f = bend.FittingParameter(5000.0,'f',1,f_prior)
w = bend.FittingParameter(10.0,'f',2,w_prior)
p = bend.FittingParameter(0.0,'f',3,p_prior)

parameters = [x0,v0,a1,b1,f,w,p]
```

Now I'll set up the model and start a chain. We have 7 parameters and very noisy data so I will start with a number of walkers on the larger side. As usual, we don't have a clue of what to expect from the chain convergence rate so it's better to use a short burn in so we can see what happens

```
model1 = bend.Model(parameters,data_t,data_x,data_x_sigma)
sampler,_,_,_ = model1.setup_sampler(500, 100, 1000)

samples, flat_samples = sampler.get_chain(), sampler.get_chain(flat=True)
label_list = model1.params_labels
bend.traceplots(samples,label_list)
bend.cornerplots(flat_samples,label_list)
```
<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb3_pic2.png" width="900">
<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb3_pic3.png" width="900">

We get a mess, but we can start cleaning things up. Notice that we get multiple peaks for the initial conditions which is terribly bad, however if we keep only the last 30% of the samples

```
bend.cornerplots(flat_samples[500*700,:],label_list)
```
<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb3_pic4.png" width="900">

We find that they are indeed approaching a single peak distribution. So we will run another chain using new priors for the initial conditions that use this result. We have 

```
x0_prior = bend.uniform_prior(-300,0)
v0_prior = bend.uniform_prior(750.0,2500.0)
x0 = bend.FittingParameter(-150.0,'x0',1,x0_prior)
v0 = bend.FittingParameter(1500.0,'v0',1,v0_prior)

parameters = [x0,v0,a1,b1,f,w,p]

model1 = bend.Model(parameters,data_t,data_x,data_x_sigma)
sampler,_,_,_ = model1.setup_sampler(500, 1000, 500)

samples, flat_samples = sampler.get_chain(), sampler.get_chain(flat=True)
label_list = model1.params_labels
bend.cornerplots(flat_samples,label_list)

solutions = [np.mean(flat_samples[:,_]) for _ in range(len(parameters))]
model1.plot_simulation(solutions)
```
<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb3_pic5.png" width="900">
<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb3_pic6.png" width="900">

Corner plots don't look good at all. And if you see the simulation with the means[^1] we can see that our damped driven oscillator is still in the initial transient state. We want our results to give an oscillator in steady state, so our chain is not going where we want it to. It's probably getting stuck in local extremas.

So we will change the strategy and will try to find a GOOD set of initial values for MCMC to begin. In order to do this we are going to try using an optimization method on our likelihood function[^2]. When we first build our MCMC model we get a function called likelihood that tells us how likely is a set of parameters for our model given the observations. When we run MCMC we try to maximize this function. But we will try a different route and see what we get. We'll begin by defining our working tools

```
from scipy.optimize import differential_evolution
target_function = model1.neg_ll
```
Differential evolution requires a set of bounds for your parameters which we'll define ad hoc from our results so far. It also gives you some optional parameters such as a seed, and two parameters concerning the foundation of this genetic algorithm the crossover (or recombination) and the mutation. We will set those as follows

```
bounds = [(-300,0),(750,2500),(-100,100),(-500,500),(1000,50000),(0,50),(-np.pi,np.pi)]
p0 = [-150,1500,2,90,5800,9.5,0.8]
rc,mut = 0.9,(1.0,1.9)
```

That set up, we can run the algorithm. I'm also using the 'workers' parameter which parallelizes the algorithm

```
diffev_solution = differential_evolution(target_function,bounds,x0=p0,mutation=mut,recombination=rc,workers=10)
```
And extract the results with

```
dev_params = diffev_solution.x
```

I got the values

```
[-1.73150188e+02,  2.27920572e+03,  1.48960440e+01,  1.81006661e+02,
         3.50990661e+04,  9.36130043e+00,  1.70851383e-01]
```

Now plotting a simulation with those parameters we get

```
model1.plot_simulation(dev_params)
```
<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb3_pic7.png" width="900">

Now this looks much more like what we are looking for. 

Back to MCMC I'll restart the chain from 0, using everything we had so far plus the results of the differential evolution algorithm. 

We have

```
# priors
x0_prior = bend.uniform_prior(-300,0)
v0_prior = bend.uniform_prior(1000.0,3000.0)
a1_prior = bend.normal_prior(0.0,1000.0)
b1_prior = bend.normal_prior(0.0,1000.0)
f_prior = bend.uniform_prior(5000.0,75000.0)
w_prior = bend.uniform_prior(0.0,30.0)
p_prior = bend.uniform_prior(-np.pi,np.pi)

# parameters
x0 = bend.FittingParameter(-173.0,'x0',1,x0_prior)
v0 = bend.FittingParameter(2280.0,'v0',1,v0_prior)
a1 = bend.FittingParameter(14.9, 'a', 1, a1_prior)
b1 = bend.FittingParameter(181.0,'b',1,b1_prior)
f = bend.FittingParameter(35100.0,'f',1,f_prior)
w = bend.FittingParameter(9.36,'f',2,w_prior)
p = bend.FittingParameter(0.171,'f',3,p_prior)

parameters = [x0,v0,a1,b1,f,w,p]

# build model
model1 = bend.Model(parameters,data_t,data_x,data_x_sigma)

# run chain

# extract and plot results

sampler,_,_,_ = model1.setup_sampler(500, 1000, 500)

samples, flat_samples = sampler.get_chain(), sampler.get_chain(flat=True)
label_list = model1.params_labels

bend.cornerplots(flat_samples,label_list)

solutions = [np.mean(flat_samples[:,_]) for _ in range(len(parameters))]
model1.plot_simulation(solutions)
```

[^1]: notice that in each corner plot the red line doesn't show the mean but the 50/50 quantile, so to get the means we calculate them from the samples using `numpy`.

[^2]: in our model we actually have the logarithm of this function, it's implemented this way so we don't get larger numbers and prevent overflows(some overflows do still happen but without this it would be worse). Also note that the function we will be optimizing is the negative log likelihood, since the optimization methods we use search for minimums. 
