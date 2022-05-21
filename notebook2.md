---
layout: custompage
permalink: /notebook2/
---

# Notebook 2: 

Remember that files summing up the work on this notebook can be found at https://github.com/tsaopy/tsaopy/tree/main/notebook2.

## Sketch of the problem

For this next problem we have the following data

<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/nb2_pic1" width="400">

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
