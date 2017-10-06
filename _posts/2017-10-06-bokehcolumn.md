---
layout: post
title: Visualising pressure swing adsorption using Bokeh
---

Recently I've been creating interactive tools for exploring data we create on the simulation of pressure swing adsorption columns to remove CO2.
These simulations were performed using [CySim][cysim], which is a adsorption cycle simulator created at The University of Edinburgh.

We simulate a four stage pressure swing adsorption cycle, operating on a gas stream which is composed of 85% nitrogen and 15% carbon dioxide,
which attempts to separate the carbon dioxide so that we can store it.
With some help from an [awesome tutorial on how to reproduce the gapminder dataset][gapminder]
I managed to produce an interactive bokeh plot of my simulation of the process.

{% include good_system.html %}

In this plot,
arrows represent gases flowing into either end of the column;
the colour of the arrows is proportional to their composition,
and the size of the arrow is proportional to the flow rate.
The lines within the column show the chemical concentration along the length of the column,
red for carbon dioxide and blue for nitrogen.
You can drag the slider to move through time in one cycle of the system,
as the system is in steady state, at the end of the cycle you have arrived back at the start of the cycle.

The cycle has four distinct stages:
In the adsorption stage, the feedstream enters the column and the carbon dioxide entering adsorbs to the solid,
allowing mostly pure nitrogen to leave the other end.
In the blowdown stage (timestep 18), the last of the nitrogen left in the column is removed by quickly pulling a vacuum on the right hand side.
A vacuum is then applied to the right hand side, which allows a mostly pure stream of carbon dioxide to leave the column.
Finally, we repressurise the column from the right hand side to prepare the column for the next cycle.

The process configuration for the above cycle (the timings of the different stages and the pressures used)
were arrived at through a genetic algorithm which optimises the
purity, recovery, energy cost and productivity of the cycle.
The details of the optimisation are probably beyond this post, but it is interesting to also visualise a "bad cycle", which I've done below.
This cycle uses the same material to adsorb CO2,
but due to being poorly set up, only achieves 72% purity and 47% recovery of the CO2, compared to 95% purity and 90% recovery above.
Using a comparison of the two visualisations, you can see why this cycle is worse than the optimised version.

{% include bad_system.html %}


[cysim]: http://www.carboncapture.eng.ed.ac.uk/lab/cysim
[bokeh]: https://bokeh.pydata.org/en/latest/
[gapminder]: https://rebeccabilbro.github.io/interactive-viz-bokeh/
