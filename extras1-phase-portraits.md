---
layout: sidepage
permalink: /phase-portraits/
---

# Extras 1: Plotting phase portraits

In this notebook we show different orders of `tsaopy` approximation for the equation

$$ \ddot{x} - \sin{(x)}\cos{(x)} + \frac{1}{2}\sin{(x)} = 0  $$

We take the original ODE and then approximate it using two different schemes (one with potential coefficients up to order 3, and another one up to order 5), and plot phase potraits using the same initial conditions, for all the three models. Here's what we got

<img src="https://raw.githubusercontent.com/tsaopy/tsaopy.github.io/main/assets/ex1/pic.png" width="700">

Notice how in this region of the phase space the 3rd order `tsaopy` model gives a rough approximation of the portrait of the original ODE, however in the 5th order model it looks almost identical.
