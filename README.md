# Tumor/T cell Microenvironment ABM

Introduction to project

## TODO
* finish README
* add references file, and somehow link these references in the processes' `defaults`
* give neighbors process units in the bounds, to be consistent with fields process.
* tumor process -- address TODOs (accessible external volume more general from timestep/diameter)
* look over Jupyter notebooks. Convert to Colab before submission -- this requires doing the pip install in the notebook.


## Vivarium

This is a vivarium project.
Visit [the Vivarium Core
documentation](https://vivarium-core.readthedocs.io/) to learn how to
use the core Vivarium engine to create computational biology models.
Check out the
[getting started](https://vivarium-core.readthedocs.io/en/latest/getting_started.html)
guide of the documentation. 


## Installation

Explain how users can get up and running.

```
pip install tumor-tcell
```

## Processes

The model is composed of two major cell types (T cells & Tumor cells), each with two separate phenotpyes. 
Each cell has an associated process. The process for each cell type contains the fundamental rules that 
govern its behavior in interacting with the other cell types and with the inputs it receives from the environment.

Testing the processes individually enables understanding if underlying parameters derived from literature values 
or primary data accurately represent behavior expected based on such research. 

### Tumor
The tumor process is focused on two states of a tumor: proliferative with low levels of immune molecules 
(MHCI and PDL1) and quiescent with high levels of immune molecules (MHCI and PDL1). Its transition from the 
proliferative state is dependent on the level of interferon gamma it is exposed to coming from the T cells. 
Both tumor types can be killed by recieving cytotoxic packets from the T cells.
<img src="/jupyter_notebooks/images/2_Tumor_process.png" alt="tumor_process" width="500"/>

### T cell
The T cell process is focused on two states of a T cell: PD1- with increased secretion of immune molecules 
(IFNg and cytotoxic packets) and PD1+ with decreased secretion of immune molecules (IFNg and cytotoxic packets). 
These immune molecules have impact of the state and death of tumor cells. Its transition from the PD1- state is 
dependent on the length of time it is engaged with tumor cells. 
<img src="/jupyter_notebooks/images/1_Tcell_process.png" alt="tcell_process" width="500"/>

## Composites

Composites are integrated models with multiple initialized processes, and whose inter-connections are specified 
by a topology. The T-cell and Tumor *processes* shown individually above are here combined with additional processes 
to create T-cell and Tumor *agents*. These include a division process, which waits for the division flag and then 
carries out division; a death process, which waits for a death flag and then removes the agent; and a local field, 
which interfaces the external environment to support uptake and secretion for each individual agent.

### Tumor Microenvironment

`TumorMicroEnvironment` is a composite model that simulates a 2D environment with agents that can move around in space, 
exchange molecules with their neighbors, and exchange molecules with a molecular field. This composite includes a 
`neighbors` process, which tracks the locations of individual agents and handles their exchanges, and a `diffusion` 
process, which operates on the molecular fields.


## Experiments
In our Multi-scale agent based model, the T cells can interact with tumor cells in the following ways:
* T cell receptor (TCR on T cells) and Major histocompatibility complex I receptor (MHCI on tumor cells) 
for activation of T cells, induction of IFNg and cytotoxic packet secretion, and slowing of T cell migration
* PD1 receptor (on T cells) and PDL1 receptor (on tumor cells) that can inhibit cell activation and induce apoptosis
* T cells secrete IFNg which tumor cells uptake and causes state switch to upregulate MHCI, PDL1, and decrease proliferation
<img src="/jupyter_notebooks/images/5_ABM.png" alt="tumor_tcell_experiment" width="500"/>


### Running the model

```
python tumor_tcell/experiments/main.py
```

