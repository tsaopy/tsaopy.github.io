---
layout: custompage
permalink: /initvals-optimization/
---

# Notebook 4: Initial values optimization

On this notebook I want to show an extra feature that `TSAOpy` has, which is finding good initial values for the MCMC chain by using an exteral optimizer. This will work best with models that you know are correct[^1] (or good enough) and you want to quickly fit the parameters to the data. To do the demonstration I simulated the following ODE which I found very interesting

[^1]: here's the problem, optimizers aren't sentient enough to be able to tell which parameters in your model are relevant or not (although MCMC is a great tool for you to figure that out), so if you just drop a bunch of parameters that may not be relevent, the optimizer may easily fall for problems like getting stuck due to correlations, getting stuck at local minimums, give you a solution that doesn't work anymore once you drop unnecessary terms, etc. 

$$ \ddot{x} -k\sin{(x)}\cos{(x)} + g \sin{(x)} = 0 $$

I'm using $k=1$, $g=0.5$, $x_0=1$, and $v_0=0.5$. The resulting data looks like this

<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb4_pic1.png" width="600">

Now, back to the ODE, I'll use the following polynomial expansion of the terms above

$$ \ddot{x} \sin{(x)}\left(g-k\cos{(x)}\right) = 0 $$

$$ \ddot{x} + (g-k)x + \frac{4k-g}{6} x^3 = 0 $$

So for our `TSAOpy` model, we will have

$$ \ddot{x} + b_1 x + b_3 x^3 = 0 $$

Now let's see how we solve this with the usual tools. I'll set up the model as

```
# priors

x0_prior = bend.uniform_prior(0.7,1.3)
v0_prior = bend.uniform_prior(0.3,0.7)
b1_prior = bend.normal_prior(0.0,10.0)
b3_prior = bend.normal_prior(0.0,10.0)
    
# parameters

x0 = bend.FittingParameter(1.0,'x0',1,x0_prior)
v0 = bend.FittingParameter(0.5,'v0',1,v0_prior)
b1 = bend.FittingParameter(0.0,'b',1,b1_prior)
b3 = bend.FittingParameter(0.0,'b',3,b3_prior)

parameters = [x0,v0,b1,b3]

# model 2

model2 = bend.VelocityModel(parameters,data_t,data_x,data_v,
                            data_x_sigma,data_v_sigma)
```
