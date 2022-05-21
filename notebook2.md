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
