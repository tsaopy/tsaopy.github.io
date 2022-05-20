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

Now for the next step I'm running the simulation. The model we will be simulating is

$$ \ddot{x} + a_1\dot{x} + b_1x = 0 $$

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
