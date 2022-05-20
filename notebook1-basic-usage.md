---
layout: custompage
permalink: /basic-usage/
---

# Notebook 1: Basic Usage

Remember that files summing up the work on this notebook can be found at https://github.com/tsaopy/tsaopy/tree/main/notebook1.

## Making some data

In this first notebook we will define a simple model, make a simulation with it, add some noise to the data, and see how to recover the original model using `TSAOpy`.


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
<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb1_pic1.png" width="400">

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
<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb1_pic2.png" width="400">

Now we are set up and we can save the data for testing `TSAOpy`. I'll do it with

```
# save results
np.savetxt('experiment_data.txt', [[t[_], x[_], v[_]] for _ in range(n_out)])
```

## Fitting the model with `TSAOpy`

Now let's forget about everything we did above. Imagine you just finished running some experiment, and you ended up with some time series that now you want to analyze. 

The first thing we need to do is to load the data. Ideally your data will be a `numpy` array, however if you want to push your luck native lists and tuples should work. To load the data I'm using 

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
<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb1_pic3.png" width="400">

Now we need to know the uncertainty of our measurements. You have two options, the simplest one is just use one value for the entire set of measurements. The other option is, if you have some way in your lab, have a unique uncertainty for each point. So you will define an uncertainty variable that is either a number representative of the uncertainty for each measurement or an array with the exact uncertainty for each variable. In this case we will just use a rough estimate for a global uncertainty and just define `data_x_sigma = 0.3`. What I tipically do is inspecting the plot, go to a zone where there are a lot of points gathered, and see how much the value of $x$ varies in a small $\Delta t$ interval.

Before we move ahead we have to decide which will be our model. Inspecting the measurement plots one sees that we have something that looks like a sinusoidal wave whose amplitud decays over time. We can also see from the graph that $x_0\approx 1$, and that since we are near an amplitude maximum near the start, then $v_0\approx 0$. Naturally in this case we will be proposing

$$ \ddot{x} + a_1\dot{x} + b_1x = 0\quad; \qquad x_0=1 \qquad v_0=0 $$ 

We ended up with the original model that we used for the simulation, and some readers may say "well it's obviously going to work if you know the right model beforehand and choose that one", and while this is true, it is also true that we didn't pick this model because we knew the solution beforehand, but because it is the simplest model that one can infer from the plot. If you are still unconvinced by this, hang on, we'll see some more examples in the following notebooks that may change your mind. That being said, let's move on.

The next thing we need is to define our parameters, for that we will be defining some parameter objects implemented speciffically for this. `TSAOpy` works with two parameter classes `FixedParameter` and `FittingParameter`. Fixed parameters will have a value (which is assumed as correct and not subject to fitting), a type, and an index. I'll talk more about those last two later on. Fitting parameters will have those same attributes, and will also have another attribute called prior, which is the probability distribution that represents our prior knowledge about that parameter. 

Defining parameters will be something like this

```
import backend as bend
p_fixed = bend.FixedParameter(1.0,'a',1)
p_variable = bend.FittingParameter(1.0,'a',1,p_prior)
```
When calling the parameter classes, the arguments are like this

1. The value. In fixed parameters this is the value we assume correct and will ALWAYS be used in simulations. In fitting parameters it will be the initial value and will be updated as the MCMC chain runs. We will tipically set the value for fitting parameters as 0, unless it's $x_0$ or $v_0$ for which we may use the first values of the time series. Another exception will be $b_1$, which corresponds to the usual harmonic potential, the linear potential term. We usually assume it's a positive number so we may set it up as 1, or some other positive value that one may infer from the plot. 
2. The type. The next argument is the parameter type which will always be a string. It's value is 'x0' and 'v0' for $x_0$ and $v_0$ respectively, 'a' for the damping terms, 'b' for the potential terms, 'c' for the mixed terms, and 'f' for the external force parameters. Don't mess these up or you will get tons of errors. 
3. The index. Indexes will be integers assigned as follows:
    1. It's always 1 for 'x0' and 'v0' parameters.
    2. For parameters 'a' and 'b' it will be the order of the term. Example, for the term $b_1x$ it will be 1, and for the term $a_3\dot{x}^3$ it will be 3.
    3. For parameters 'c' it will be a pair of indices in the form of a touple. The first value will be the order of the position factor, and the second value will be the order of the velocity factor. Example, if the term is $c_{21}x^2\dot{c}$ the indices will be (2,1), and if the term is $c_{23}x^2\dot{c}^3$ the indices will be (2,3).
    4. Finally for the external or driving force parameters we use 1 for $F_0$, 2 for $\omega$, and 3 for $\phi$. 
