# jInv Examples

The folder `2016-RTH-SISC` contains code examples used to generate the results in the [paper](http://arxiv.org/abs/1606.07399).  

## Inversion Examples

For educational purposes we have created some examples and inversion tests. These are typically unrealistically small-scale but show how jInv can be used and set up to solve also bigger problems. The examples are also a great starting point for experiments with different parameters, solvers, etc.  Currently, there are the following examples:

1. `exDCResistivity.ipynb` - 3D example for parameter estimation for the Poisson equation. The example is motivated by the geophysical imaging technique of DC Resistivity.
1. `exEikonal.ipynb` - 3D example for parameter estimation for the nonlinear Eikonal equation. Example is motivated by travel time tomography, which is a geophysical imaging technique.
1. `exJointEikonalDC.ipynb` - example of a multiphysics inversion. The example combines both previous datasets and aims at estimating a single model that explains both measurements.  
