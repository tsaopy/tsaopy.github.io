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

We find that they are indeed approaching single peak distributions but there is a lot of noise and outliers. So I will start another chain with some new priors and initial values using this result. I will not update $a_1$, $b_1$, and $F_0$ priors at this point because I want to allow those parameters to run all over the place. We have 

```
x0_prior = bend.uniform_prior(-300,-20)
v0_prior = bend.uniform_prior(750,3000.0)
a1_prior = bend.normal_prior(0.0,1000.0)
b1_prior = bend.normal_prior(0.0,1000.0)
f_prior = bend.uniform_prior(0.0,50000.0)
w_prior = bend.normal_prior(9.5,1.0)
p_prior = bend.normal_prior(0.75,0.3)

x0 = bend.FittingParameter(-160.0,'x0',1,x0_prior)
v0 = bend.FittingParameter(1550.0,'v0',1,v0_prior)
a1 = bend.FittingParameter(2.0, 'a', 1, a1_prior)
b1 = bend.FittingParameter(90.0,'b',1,b1_prior)
f = bend.FittingParameter(5900.0,'f',1,f_prior)
w = bend.FittingParameter(9.4,'f',2,w_prior)
p = bend.FittingParameter(0.75,'f',3,p_prior)

parameters = [x0,v0,a1,b1,f,w,p]

# model
model1 = bend.Model(parameters,data_t,data_x,data_x_sigma)
sampler,_,_,_ = model1.setup_sampler(500, 1000, 500)

samples, flat_samples = sampler.get_chain(), sampler.get_chain(flat=True)
label_list = model1.params_labels

bend.cornerplots(flat_samples,label_list)
```
<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb3_pic5.png" width="900">

Now we are getting considerably better posteriors. However, notice how there's a tail in the $a_1$ vs $F_0$ corner plot, and there's still some outliers scattered on the other plots. We want to get rid of that and focus on the good looking peaks we got, and for that we will try running yet another chain with updated priors. 


```
x0_prior = bend.normal_prior(-148,31)
v0_prior = bend.normal_prior(1920,380.0)
a1_prior = bend.normal_prior(6.0,5.0)
b1_prior = bend.normal_prior(90.0,20.0)
f_prior = bend.normal_prior(12000.0,1000.0)
w_prior = bend.normal_prior(9.40,0.01)
p_prior = bend.normal_prior(0.63,0.31)

x0 = bend.FittingParameter(-148.0,'x0',1,x0_prior)
v0 = bend.FittingParameter(1920.0,'v0',1,v0_prior)
a1 = bend.FittingParameter(6.0, 'a', 1, a1_prior)
b1 = bend.FittingParameter(90.0,'b',1,b1_prior)
f = bend.FittingParameter(12000.0,'f',1,f_prior)
w = bend.FittingParameter(9.40,'f',2,w_prior)
p = bend.FittingParameter(0.63,'f',3,p_prior)

parameters = [x0,v0,a1,b1,f,w,p]

# model
model1 = bend.Model(parameters,data_t,data_x,data_x_sigma)
sampler,_,_,_ = model1.setup_sampler(500, 1000, 500)

samples, flat_samples = sampler.get_chain(), sampler.get_chain(flat=True)
label_list = model1.params_labels

bend.cornerplots(flat_samples,label_list)
```
