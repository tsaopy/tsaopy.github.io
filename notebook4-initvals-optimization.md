---
layout: sidepage
permalink: /initvals-optimization/
---

Remember that files summing up the work on this notebook can be found at https://github.com/tsaopy/tsaopy/tree/main/notebook2.

# Notebook 4: Initial values optimization

On this notebook I want to show an extra feature that `TSAOpy` has, which is finding good initial values for the MCMC chain by using an exteral optimizer. This will work best with models that you know are correct[^1] (or good enough) and you want to quickly fit the parameters to the data. To do the demonstration I simulated the following ODE which I found very interesting

[^1]: here's the problem, optimizers aren't sentient enough to be able to tell which parameters in your model are relevant or not (although MCMC is a great tool for you to figure that out), so if you just drop a bunch of parameters that may not be relevent, the optimizer may easily fall for problems like getting stuck due to correlations, getting stuck at local minimums, give you a solution that doesn't work anymore once you drop unnecessary terms, etc. 

$$ \ddot{x} -k \sin{(x)}\cos{(x)} + g \sin{(x)} = 0 $$

I'm using $k=1$, $g=0.5$, $x_0=1$, and $v_0=0.5$. The resulting data looks like this

<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb4_pic1.png" width="600">

Now, back to the ODE, I'll use the following polynomial expansion of the terms above

$$ \ddot{x} + \sin{(x)}\left(g-k \cos{(x)}\right) = 0 $$

$$ \ddot{x} + (g-k)x + \frac{4k-g}{6} x^3 = 0 $$

So for our `TSAOpy` model, we will have

$$ \ddot{x} + b_1 x + b_3 x^3 = 0 $$

I'll set up the code model as

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

Now let's see how we solve this with the usual procedure. We'd tipically start by running a preliminary chain to see where things go from zero

```
sampler,_,_,_ = model2.setup_sampler(300, 50, 500)
samples, flat_samples = sampler.get_chain(), sampler.get_chain(flat=True)

label_list = model2.params_labels
bend.cornerplots(flat_samples,label_list)
bend.traceplots(samples,label_list)
```
<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb4_pic2.png" width="600">
<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb4_pic3.png" width="600">

And then say "well it's taking roughly 300 steps to start converging to some peaks which now we know their value"[^2]

[^2]: remember that the easiest way to calculate those peaks is dropping the samples from the first 300 steps and then taking the means of the remaning samples, like this `peaks = [np.mean(flat_samples[300*300:,_]) for _ in range(len(parameters))]`.

So we update priors and start a new chain like this 

```
# priors
x0_prior = bend.normal_prior(1.0,0.1)
v0_prior = bend.normal_prior(0.44,0.8)
b1_prior = bend.normal_prior(-0.32,0.2)
b3_prior = bend.normal_prior(0.26,0.2)
    
# parameters
x0 = bend.FittingParameter(1.0,'x0',1,x0_prior)
v0 = bend.FittingParameter(0.44,'v0',1,v0_prior)
b1 = bend.FittingParameter(-0.32,'b',1,b1_prior)
b3 = bend.FittingParameter(0.26,'b',3,b3_prior)

parameters = [x0,v0,b1,b3]

# model 2

model2 = bend.VelocityModel(parameters,data_t,data_x,data_v,
                            data_x_sigma,data_v_sigma)

sampler,_,_,_ = model2.setup_sampler(300, 500, 500)
samples, flat_samples = sampler.get_chain(), sampler.get_chain(flat=True)

label_list = model2.params_labels
bend.cornerplots(flat_samples,label_list)
bend.traceplots(samples,label_list)
```
<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb4_pic4.png" width="600">
<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb4_pic5.png" width="600">

And we got the job done by running at least 3 chains of 500 steps each (you may need more chains if you don't choose adequate walker and step numbers).

Now, let's go back to square one and try the alternative procedure. We have our model

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

We want to optimize the initial values before running any chains. For that we have to optimize the negative logarithmic likelihood function[^3], which we get from the model object as

[^3]: the likelihood function is a function that, given a model and a set of parameters, tells you how likely it is that those parameters explain some data. Finding the maximum of the likelihood function gives you the "best fitting" or more likely parameters, and this is what MCMC does. However one may try to maximize this function by any other means and this is what we are about to do. However, notice that I mentioned negative logarithmic likelihood. First of all we work with the logarithm of the likelihood function because statisticians say it has some nice properties that make things work better. We accept it with no arguments and move on. Then, we work with the negative log likelihood because we need to maximize it, and most optimizers search for minimums.

```
neg_ll = model2.neg_ll
```

Now, I will use `scipy.optimize.differential_evolution` for the optimization. You can pick any other optimizer that you like but I found this one to work OK. For this method we need a set of bounds for the parameters, I'll pick some bounds to roughly mimic the priors

```
from scipy.optimize import differential_evolution
bounds = [(0.7,1.3),(0.3,0.7),(-30,30),(-30,30)]
```

Differential Evolution has many parameters, of which I'm mostly interested in mutation, population size, and max number of iterations. Mutation and population size are typical parameters of any genetic algorithm, in this case I'm picking values that, compared with the defaults, will make the optimizer more stable and consistent but slower. I will allow the algorithm to run for a larger number of iterations, since for this kind of problem the default maximum is not enough. You can also use the workers parameter to parallelize the computations.

Summing up, this is what we have

```
diffev_solution = differential_evolution(neg_ll,bounds=bounds,popsize=50,
                              mutation=(1.0,1.9),maxiter=2000,workers=6)

new_initvals = diffev_solution.x
```
We can try plotting a simulation with the new values that we just found

```
model2.plot_simulation(new_initvals)
```
<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb4_pic6.png" width="600">

The results look very good, so we are going to update the model's initial values and do the MCMC run

```
model2.update_initvals(new_initvals)
sampler,_,_,_ = model2.setup_sampler(300, 100, 500)
samples, flat_samples = sampler.get_chain(), sampler.get_chain(flat=True)
label_list = model2.params_labels
bend.cornerplots(flat_samples,label_list)
```
<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb4_pic7.png" width="600">

And we get identical posteriors than what we got before, but now we only had to run one chain with a small burn in. 
