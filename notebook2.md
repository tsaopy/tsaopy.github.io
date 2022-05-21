---
layout: custompage
permalink: /notebook2/
---

# Notebook 2: 

Remember that files summing up the work on this notebook can be found at https://github.com/tsaopy/tsaopy/tree/main/notebook2.

### Note: for some reason the page layout wouldn't let me to show bigger pictures so you may want to right click some of them and look them at full size. 

## Sketch of the problem

For this next problem we have the following data

<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb2_pic1.png" width="600">

From the plots we can't say much more than that the initial conditions are roughly $x_0\approx2$ and $v_0\approx-0.3$. I'm also going to estimate the uncertainties as

```
data_x_sigma,data_v_sigma = 0.3,0.4
```

Now we have to pick a model. There's really nothing indicative of what's going on there so I will propose a model with all parameters up to order 2. Hopefully, parameters which don't really fit this model nicely will converge to 0. The ODE will look like this

$$ \ddot{x} + a_1\dot{x} + a_2\dot{x}|\dot{x}| + b_1x + b_2x^2 + c_{11}x\dot{x} + c_{12}x\dot{x}^2 + c_{21}x^2\dot{x} + c_{22}x^2\dot{x}^2 = 0 $$

We don't know anything about the parameters besides rough estimations of the initial conditions, so I will pick normal priors like this

```
# priors
x0_prior = bend.normal_prior(2.0,5.0)
v0_prior = bend.normal_prior(-0.3,5.0)
a1_prior = bend.normal_prior(0.0,10.0)
a2_prior = bend.normal_prior(0.0,10.0)
b1_prior = bend.normal_prior(0.0,10.0)
b2_prior = bend.normal_prior(0.0,10.0)
c11_prior = bend.normal_prior(0.0,10.0)
c12_prior = bend.normal_prior(0.0,10.0)
c21_prior = bend.normal_prior(0.0,10.0)
c22_prior = bend.normal_prior(0.0,10.0)
```
Then I will define the parameter objects and the list of parameters as follows

```
# parameters
x0 = bend.FittingParameter(2.0,'x0',1,x0_prior)
v0 = bend.FittingParameter(-0.3,'v0',1,v0_prior)
a1 = bend.FittingParameter(0.0, 'a', 1, a1_prior)
a2 = bend.FittingParameter(0.0, 'a', 2, a2_prior)
b1 = bend.FittingParameter(0.0,'b',1,b1_prior)
b2 = bend.FittingParameter(0.0,'b',2,b2_prior)
c11 = bend.FittingParameter(0.0,'c',(1,1),c11_prior)
c12 = bend.FittingParameter(0.0,'c',(1,2),c12_prior)
c21 = bend.FittingParameter(0.0,'c',(2,1),c21_prior)
c22 = bend.FittingParameter(0.0,'c',(2,2),c22_prior)

parameters = [x0,v0,a1,a2,b1,b2,c11,c12,c21,c22]
```
Now we can build our model object. `TSAOpy` comes with a second type of model which also fits parameters to $v(t)$ measurements. We set it up like this

```
model2 = bend.VelocityModel(parameters,data_t,data_x,data_v,
                            data_x_sigma,data_v_sigma)
```

Now it is time to begin running chains. Since this problem is considerably more complex than the one on the previous notebook, I'll begin by starting with a higher number of walkers. Using a higher number of walkers helps the chain converge to the right place in less steps, and sometimes helps prevents some problems like some walkers not moving, or converging to different values than the rest of the walkers. 

About the chain steps, let's remember that burn in steps are thrown away. This is important because if you have some intuition about the number of steps your chain will need to converge, you can just set that length, plus some extra just in case, for your burn in phase, and then for the production phase your chain already converged so you are really just drawing samples from your correct solution. So if your burn in phase is long enough for your chain to already have converged to the target posterior by it's end, and your production phase is long enough for you to collect enough samples according to your criteria, then in just one run you will get the correct result. 

Now, what if you have no clue about how long will it take for your chain to converge? Remember that we don't even know if this a good model at all. So now we don't want to drop samples, we want to see what the walkers are doing. So we will run a very short burn in and have a longer production so we can see where the chain is going. I'm setting up the sampler as

```
sampler,_,_,_ = model2.setup_sampler(500, 50, 500)
```
And after the run is finished the results we care about are
```
samples, flat_samples = sampler.get_chain(), sampler.get_chain(flat=True)
label_list = model2.params_labels
bend.cornerplots(flat_samples,label_list)
bend.traceplots(samples,label_list)
```
<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb2_pic2.png" width="900">
<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb2_pic3.png" width="900">

Notice the following
1. In the trace plots we are seeing the dreaded problem of many walkers being stuck. Remember that we want each individual walker to be constantly changing, and a group behaviour of oscillating around a value, but we don't want individual walkers to show a straight horizontal line in the trace plot. That means things aren't going anywhere. 
2. This also can be seen in the corner plots where we have points spreaded all over the place instead of being centered on a high density zone. What's happening is that many walkers are stuck in a value of their own and not following the others.

With that in mind, we are going to try running another chain, longer, and with more walkers. So now we try this

```
sampler,_,_,_ = model2.setup_sampler(1000, 50, 1000)
samples, flat_samples = sampler.get_chain(), sampler.get_chain(flat=True)
bend.cornerplots(flat_samples,label_list)
bend.traceplots(samples,label_list)
```
<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb2_pic4.png" width="900">
<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb2_pic5.png" width="900">

We get more of the same. Let's see what happens if we plot only the last 300 steps
```
bend.cornerplots(flat_samples[1000*700:],label_list)
bend.traceplots(samples[700:],label_list)
```
<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb2_pic6.png" width="900">
<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb2_pic7.png" width="900">

We can see that many parameters have very sharp peaks around a value, but there are still many walkers that aren't converging there. So now we can try something different. We will build another model, but we will update the priors to use this results. We have

```
# priors
x0_prior = bend.normal_prior(1.9,0.2)
v0_prior = bend.normal_prior(0.3,0.4)
a1_prior = bend.normal_prior(-1.9,2.0)
a2_prior = bend.normal_prior(0.25,0.3)
b1_prior = bend.normal_prior(0.9,0.3)
b2_prior = bend.normal_prior(0.0,1.0)
c11_prior = bend.normal_prior(0.0,1.0)
c12_prior = bend.normal_prior(0.0,1.0)
c21_prior = bend.normal_prior(1.6,1.5)
c22_prior = bend.normal_prior(0.0,1.0)

# parameters
x0 = bend.FittingParameter(1.9,'x0',1,x0_prior)
v0 = bend.FittingParameter(0.3,'v0',1,v0_prior)
a1 = bend.FittingParameter(-1.9, 'a', 1, a1_prior)
a2 = bend.FittingParameter(0.25, 'a', 2, a2_prior)
b1 = bend.FittingParameter(0.9,'b',1,b1_prior)
b2 = bend.FittingParameter(0.0,'b',2,b2_prior)
c11 = bend.FittingParameter(0.0,'c',(1,1),c11_prior)
c12 = bend.FittingParameter(0.0,'c',(1,2),c12_prior)
c21 = bend.FittingParameter(1.6,'c',(2,1),c21_prior)
c22 = bend.FittingParameter(0.0,'c',(2,2),c22_prior)

parameters = [x0,v0,a1,a2,b1,b2,c11,c12,c21,c22]
```
Notice that some parameters had a posterior like $0.07\pm0.15$. A posterior like that is practically impossible to distinguish from 0, it's very likely that those parameters do not correspond to the model and may converge to zero, so I'll keep them as zero but reduce the standard deviation since they weren't spread over so thich intervals.

And now we build the model again, and run another chain

```
model2 = bend.VelocityModel(parameters,data_t,data_x,data_v,
                            data_x_sigma,data_v_sigma)

sampler,_,_,_ = model2.setup_sampler(1000, 50, 1000)
samples, flat_samples = sampler.get_chain(), sampler.get_chain(flat=True)

bend.cornerplots(flat_samples,label_list)
bend.traceplots(samples,label_list)
```
<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb2_pic8.png" width="900">
<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb2_pic9.png" width="900">

Now this looks much better. Notice that $b_2$, $c_{11}$, $c_{12}$, and $c_{22}$, are practically 0. Their posteriors all have means close to 0, and have SDs of absolute magnitude bigger than the mean itself, that's why we really can't distinguish this value from 0. So we will try running another chain, this time dropping those parameters. 

I'll reduce walkers a bit since we have simplified the model, and will use a burn in phase and see if we could get a chain converging to a good looking posterior.

```
parameters = [x0,v0,a1,a2,b1,c21]

model2 = bend.VelocityModel(parameters,data_t,data_x,data_v,
                            data_x_sigma,data_v_sigma)

sampler,_,_,_ = model2.setup_sampler(500, 1000, 500)
samples, flat_samples = sampler.get_chain(), sampler.get_chain(flat=True)

bend.cornerplots(flat_samples,label_list)
```
<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb2_pic10.png" width="900">

Now we got a set of posteriors that look from a chain that has indeed converged. The next thing to notice is how badly $a_2$ is correlated to $a_1$ and $b_1$. We also point out that it's mean value went from roughly 0.26 to 0.05, and now it's very similar to 0. Because of this I also want to drop $a_2$ and see what we get. 

```
parameters = [x0,v0,a1,b1,c21]

model2 = bend.VelocityModel(parameters,data_t,data_x,data_v,
                            data_x_sigma,data_v_sigma)

sampler,_,_,_ = model2.setup_sampler(500, 1000, 500)
samples, flat_samples = sampler.get_chain(), sampler.get_chain(flat=True)

bend.cornerplots(flat_samples,label_list)
```
<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb2_pic11.png" width="900">

Unexpectedly I found some noisy bimodal sort of behavior for the $x_0$ and $v_0$ posteriors, so we try again after updating their priors

```
parameters = [x0,v0,a1,b1,c21]

model2 = bend.VelocityModel(parameters,data_t,data_x,data_v,
                            data_x_sigma,data_v_sigma)

sampler,_,_,_ = model2.setup_sampler(500, 1000, 500)
samples, flat_samples = sampler.get_chain(), sampler.get_chain(flat=True)

bend.cornerplots(flat_samples,label_list)
```
