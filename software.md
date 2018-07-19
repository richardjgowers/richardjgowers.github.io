---
layout: page
title: Software
permalink: /software/
---

I enjoy writing software (mostly Python, some C/Fortran, some other things)
to solve the problems in chemistry and engineering
for my research.  Here's a semi complete list of things I have worked on.

## Software I make

### [MDAnalysis](http://mdanalysis.org)
> Analyse all your MD results in Python!

MDAnalysis is a Python package which lets you read and write
molecular simulation data.  I've been working on MDAnalysis since
2013(!) and it is how I fell down the rabbit hole of open source software.

### [GCMCWorkflow](https://github.com/richardjgowers/GCMCworkflow)
> Automate your material screening!

GCMCWorkflow is a workflow tool (built on top of Fireworks) for running Monte Carlo simulations
of gas adsorption.  It doesn't run the simulations (it wraps around 
[Raspa](https://www.iraspa.org/RASPA/index.html)), but instead lets
you write large complex workflows to do complex tasks.

### [Datreant](http://datreant.org)
> Filesystem as a database!

Datreant lets you view your file system as a database, allowing you
to filter, sort and search through your data.

### [Hydraspa](https://github.com/richardjgowers/hydraspa)

A bunch of utilities for Raspa and the CoreMof database.
Create/manipulate system inputs and parse outputs.

### [pmda](https://www.mdanalysis.org/pmda/)

Allows MDAnalysis to be used in massively parallel ways.

### [MDAPackmol](https://github.com/MDAnalysis/MDAPackmol)

A wrapper around [Packmol](http://m3g.iqm.unicamp.br/packmol/home.shtml)
which lets you pass in MDAnalysis objects into Packmol.
This then makes system creation for a variety of simulation programs
easier as topology information (bonds etc) is correctly preserved.

### [IBIsCO](https://github.com/richardjgowers/ibisco)

A molecular simulation program I developed during my PhD.
It is capable of running dual scale/hybrid scale simulations, ie
both atomistic and coarse-grained particles at the same time.

## Software I've contributed to

> (Stuff which isn't really mine, but I liked and used enough to fix things)

- Fireworks
- Inspyred
