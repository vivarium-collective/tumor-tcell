# Tumor/T cell Microenvironment ABM

Introduction to project


## Installation

Explain how users can get up and running.

```
pip install tumor-tcell
```

## Processes

Briefly describe these.
* Tcell
* Tumor
* Neighbors

## Composites

Composites are integrated models with multiple initialized processes, and whose inter-connections 
are specified by a topology. The T-cell and Tumor *agents* include a division 
process, which waits for the division flag and then carries out division; a death process, which 
waits for a death flag and then removes the agent; and a local field, which interfaces the external 
environment to support uptake and secretion for each individual agent.

## Running the model

```
python tumor_tcell/experiments/main.py
```

## Vivarium

This is a vivarium project.
Visit [the Vivarium Core
documentation](https://vivarium-core.readthedocs.io/) to learn how to
use the core Vivarium engine to create computational biology models.
Check out the
[getting started](https://vivarium-core.readthedocs.io/en/latest/getting_started.html)
guide of the documentation. 

------------------------------------------------------------------------
