---
layout: sidepage
permalink: /phase-portraits/
---

# Extras 1: Plotting phase portraits

Following [notebook 4](https://tsaopy.github.io/initvals-optimization/), I thought it would be interesting to make some phase portraits to compare the performance of different models optimized with `tsaopy`.

So I picked the original ODE and then approximated it using two different schemes (one of them is the model used on the notebook), and plotted phase potraits using the same initial conditions, for all the three models. Here's what I got

<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/ex1/pic.png" width="700">

Notice that in this region of the phase space the 3rd order `tsaopy` model gives a rough approximation of the portrait of the original ODE, however in the 5th order model it looks almost identical.
