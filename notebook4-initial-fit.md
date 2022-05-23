---
layout: custompage
permalink: /initial-fit/
---

# Notebook 4: Bead on a hoop

On this notebook I want to show an extra feature that `TSAOpy` has, which is finding good initial values for the MCMC chain by using an exteral optimizer. And to use as demonstration I will simulate the following ODE which I found very interesting for reasons I'll be explaining shortly

$$ \ddot{x} -k\sin{(x)}\cos{(x)} + g \sin{(x)} = 0 $$
