---
layout: sidepage
permalink: /basic-usage/
---

# Notebook 1: Basic Usage

{:toc}

## Making some data

In this first notebook we will define a simple model, make a simulation with it, add some noise to the data, and see how to recover the original model using `tsaopy`.


I'll start by setting up the numerical integrator (I'm gonna skip the details since I'm assuming you are familiar with writing a simple Runge-Kutta integrator). I implemented this as

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

Now for the next step we are running the simulation. The model we will be simulating is

$$ \ddot{x} + a_1\dot{x} + b_1x = 0\quad; \qquad x_0=1 \qquad v_0=0 $$

```
# simulation
a1,b1 = 0.3,1.0
deriv = lambda X : np.array([  X[1],  -a1*X[1]-b1*X[0]  ])

X0 = np.array([1.0,0.0])

result = solve_ivp(deriv,X0,(0,30),0.01)

t,x,v = result[:,0],result[:,1],result[:,2]
```

We get as result 3 arrays with the t, x, and v data. Then I'm going run the "fix time series length" blocks in order to set a specific number of points for the time series. 

```
# fix time series length
n_data = len(t)
n_out = 1000
step = n_data//n_out

x, v, t = x[::step], v[::step], t[::step]
while len(t) > n_out:
    x, v, t = x[:-1], v[:-1], t[:-1]
```

Now let's take a moment to plot the results and see what we have so far


```
# plot
plt.figure(figsize=(7, 5), dpi=150)
plt.plot(t, x, color = 'tab:red', label='position')
plt.plot(t, v, color='tab:purple', label='velocity')
plt.legend()
plt.show()
```
<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb1/nb1_pic1.png" width="700">

Now we will add some noise to the data and plot the results.

```
# aux noise function
np.random.seed(205)
u_noise,n_noise = 1e-1,1e-1
noise = lambda : np.random.uniform(-u_noise,u_noise) +\
    np.random.normal(0,n_noise)

# add noise
for i in range(n_out):
    x[i] = x[i] + noise()*0.3
    v[i] = v[i] + noise()*0.2
    
# plot
plt.figure(figsize=(7, 5), dpi=150)
plt.plot(t, x, color = 'tab:red', label='position')
plt.plot(t, v, color='tab:purple', label='velocity')
plt.legend()
plt.show()
```
<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb1/nb1_pic2.png" width="700">

Now we are set up and we can save the data for testing `tsaopy`. I'll do it with

```
# save results
np.savetxt('experiment_data.txt', [[t[_], x[_], v[_]] for _ in range(n_out)])
```

## Fitting the model with TSAOpy

Now let's forget about everything we did above. Imagine you just finished running some experiment, and you ended up with some time series that now you want to analyze. 

The first thing we need to do is to load the data. Ideally your data will be a `numpy` array, however native lists and tuples should work. To load the data I'm using 

```
import numpy as np

# load data
data = np.loadtxt('experiment_data.txt')
data_t,data_x = data[:,0],data[:,1]
```

Let's make a plot to see what we have to work with. 
```
# plot
plt.figure(figsize=(7, 5), dpi=150)
plt.scatter(data_t, data_x, color = 'tab:red', s=0.5, label='x(t)')
plt.legend()
plt.show()
```

<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb1/nb1_pic3.png" width="700">

Now we need to know the uncertainty of our $x$ measurements. You have two options, the simplest one is just use one value for the entire set of measurements. The other option is, if you have some way to do it in your lab, save a unique value of uncertainty for each point. So you will define an uncertainty variable that is either a number representative of the uncertainty for all measurements or an array with a unique value for each point. In this case we will just use a rough estimate for all points and just define `data_x_sigma = 0.15`.

Before we move ahead we have to decide which will be our model. Inspecting the measurement plots one sees that we have something that looks like a sinusoidal wave whose amplitud decays over time. We can also see from the graph that $x_0\approx 1$, and that since we are near an amplitude maximum near the start, then $v_0\approx 0$. Naturally in this case we will be proposing

$$ \ddot{x} + a_1\dot{x} + b_1x = 0\quad; \qquad x_0=1 \qquad v_0=0 $$ 

The next thing we need is to define our `tsaopy` parameter objects. `tsaopy` works with two parameter classes `Fixed` and `Fitting`. Fixed parameters will have a value (which is assumed as correct and not subject to fitting), a ptype, and an index. I'll talk more about those last two later on. Fitting parameters will have those same attributes, and will also have another attribute called prior, which is the probability distribution that represents our prior knowledge about that parameter. 

Defining parameters will be something like this

```
import tsaopy
p_fixed = tsaopy.parameters.Fixed(1.0,'a',1)
p_variable = tsaopy.parameters.Fitting(1.0,'a',1,p_prior)
```
When calling the parameter classes, these are the arguments

1. The value. In fixed parameters this is the value we assume correct and will ALWAYS be used in simulations. In fitting parameters it will be the initial value and will be updated as the MCMC chain runs. 
2. The type. The next argument is the parameter ptype which will always be a string. It's value is 'x0' and 'v0' for $x_0$ and $v_0$ respectively, 'a' for the damping terms, 'b' for the potential terms, 'c' for the mixed terms, and 'f' for the external force parameters. Don't mess these up or you will get tons of errors. 
3. The index. Index is set as None as default, but must be set for 'a', 'b', 'c', and 'f' parameters or you will get errors. The general rules are:
    1. No need to change it from the default None for 'ep', 'x0', and 'v0' parameters.
    2. For parameters 'a' and 'b' it will be the order of the term. Example, for the term $b_1x$ it will be 1, and for the term $a_3\dot{x}^3$ it will be 3.
    3. For parameters 'c' it will be a pair of indices in the form of a touple. The first value will be the order of the position factor, and the second value will be the order of the velocity factor. Example, if the term is $c_{21}x^2\dot{x}$ the indices will be (2,1), and if the term is $c_{13}x\dot{x}^3$ the indices will be (1,3).
    4. Finally for the external or driving force parameters we use 1 for $F_0$, 2 for $\omega$, and 3 for $\phi$. 

A somewhat obvious remark at this point should be that in all models we should have at least 2 parameters corresponding to the initial conditions, and at least one term in the ODE, otherwise we will always get straight lines as a result. 

Now, back at the problem, our model has 4 parameters, $x_0$, $v_0$, $a_1$, and $v_1$. We will not fix any of them so we need a prior for each of them. We have some prior classes in the backend, in this case I'll use the uniform prior, which takes the min and max values as argument to create the object. The code is going to be

```
# priors
x0_prior = tsaopy.tools.uniform_prior(0.7,1.3)
v0_prior = tsaopy.tools.uniform_prior(-1.0,1.0)
a1_prior = tsaopy.tools.uniform_prior(-5.0,5.0)
b1_prior = tsaopy.tools.uniform_prior(0.0,5.0)
```
We know that $x_0$ and $v_0$ are roughly 1 and 0 respectively, so we just take an interval centered at the value we think it's correct. For $a_1$ we know nothing so we just make a larger interval centered at 0. And finally for $b_1$ we know that it is possitive so we take an interval starting from 0. Now we can define out parameter objects and we have

```
# parameters
x0 = tsaopy.parameters.Fitting(1.0,'x0',1,x0_prior)
v0 = tsaopy.parameters.Fitting(0.0,'v0',1,v0_prior)
a1 = tsaopy.parameters.Fitting(0.0, 'a', 1, a1_prior)
b1 = tsaopy.parameters.Fitting(0.5,'b',1,b1_prior)

parameters = [x0,v0,a1,b1]
```
notice that I also saved them on a list. So now we have our model, our data, and our priors, and we wraped it as required by `tsaopy`. The next step is building the `tsaopy` model object which will condense everything in a single object. We call it with 

```
model1 = tsaopy.models.PModel(parameters,data_t,data_x,data_x_sigma)
```

This object will have some QOL methods such as `model.plot_measurements()` to plot your data and check that it was correctly loaded, `model.update_initvals(p0)` which allows you to enter a new set of initial values in case you want to, and others which will probably be described in a future API. During the development of the notebooks I will be mentioning the most useful ones anyways. 

Now once we defined the `tsaopy` Model, we will use an intrinsic method of this class to build an `emcee` Sampler to use the MCMC method, and start running it. The line will be something like

```
sampler,_,_,_ = model1.setup_sampler(200, 300, 300)
```

The first argument is the number of walkers, the second one the number of burn in steps, and the third one the number of production steps. The number of walkers is tipically set at few hundreds, the more you use the faster the chain may converge, but it will also take more function evaluations per step and will take longer per step. The number of burn in steps is the ammount of steps that you think the chain will take to converge. Running this line will also get MCMC running so you are in for a little wait until it's finished. After the `emcee` run is finished we will call the next line which will extract the sample chains from the sampler. We will run

```
samples, flat_samples = sampler.get_chain(), sampler.get_chain(flat=True)
label_list = model1.params_labels
```

The only difference between samples and flat samples is that samples stores each step of the chain separatedly for each walker, and flat samples will just throw all samples together regardless of which walker it came from. The diffence is that some of the plots that we will make know take different sample formats as arguments. The labels call is just to store the name of each parameter. 

Now we will make three plots to show the results. The first one will be a corner plot which goes as follows

```
tsaopy.tools.cornerplots(flat_samples,label_list)
```

<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb1/nb1_pic4.png" width="700">

Here we get the posterior distribution for each parameter, and plots of the posteriors for each pair of parameters showing possible correlations.

The next two plots will be useful to analyze wether the method has converged or not. The following plot is called the trace plot and will show the value of each parameter for each walker at each time step. 

```
tsaopy.tools.traceplots(samples,label_list)
```

<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb1/nb1_pic5.png" width="700">

A first diagnose the trace plot gives is that if for any of the parameters the 'cloud' is not exacly horizontal and uniform but it's bent, or shows thinning or thickening, then when can't say for sure the chain converged. A first attempt at fixing this is running another chain with a longer burn in phase. Another problem we may encounter is finding single perfectly horizontal lines (as in individual lines, not the cloud). This means that this walker did not change it's value at all, perhaps because it got stuck in a local extrema, or because we didn't define priors properly. We will talk more about these problems later on. 

The next plot is the autocorrelation plot. It's very simple. If all of the walkers for a parameter converged then the function quickly drops from 1 to 0 and remains oscillating around 0. If it doesn't drop to 0 quicky, let's say before the first third of the chain, then we can't say for sure it converged. Again the first attempt to fix this should be running a longer chain (this time take longer burn in and production phases). 

It goes like this
```
tsaopy.tools.autocplots(flat_samples,label_list)
```

<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb1/nb1_pic6.png" width="700">

Now given what we just said about the trace and autocorrelation plots, and the plots we got, we can assume that the chain has indeed converged. 

Finally we will plot a simulation using as values the mean of the posterior for each parameter, and compare it to the data we had. 

```
solutions = [np.mean(flat_samples[:,_]) for _ in range(len(parameters))]
model1.plot_simulation(solutions)
```

<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb1/nb1_pic7.png" width="700">

### Practice ideas

1. Run another chain with only one step of burn in and check how the corner plots and trace plots look. This will show you how the parameter space is explored during the burn in.
2. Run more chains with different values for the data uncertainty and check what you get.
3. Run more chains with different sets of fixed initial conditions and check what you get. 
4. Try reducing the length of the burn in phase to find the minimum length for which the chain has converged.
5. Run a chain, and use the results as priors and initial values for a new chain (you should use the tsaopy.tools.normal_prior(x0,sigma) object for the new priors). Check that this will allow you to run a shorter burn in phase in the new chain.
6. Run a chain with more free parameters and check the corner plots for correlations. Is the chain converging as easily as before? What value are the free parameters converging to? Note: I added the $a_2$ and $b_2$ terms and with 200 walkers now it takes roughly 500 steps for the chain to converge. 
