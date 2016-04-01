[![Build Status](https://travis-ci.org/lruthotto/jInv.jl.svg?branch=master)](https://travis-ci.org/lruthotto/jInv.jl) 
[![Coverage Status](https://coveralls.io/repos/github/lruthotto/jInv.jl/badge.svg?branch=master)](https://coveralls.io/github/lruthotto/jInv.jl?branch=master)
[![Build status](https://ci.appveyor.com/api/projects/status/0pxgtmm08b0w6wgh?svg=true)](https://ci.appveyor.com/project/lruthotto/jinv-jl-81lel)

*This fork is used to restructure some things in jInv.InverseSolve. The goal is to make the package easier to use for multiphysics problems and random sampling of data. Making this easier, required some major changes in the way InverseParam is used. Once everything is stable, this fork will be merged. Please look at the difference between this and the original version and let me know if you have any comments.*

# jInv

`jInv` is a flexible framework for PDE parameter estimation in Julia. It provides easy to extend core functions used in PDE-constrained inverse problems.

# Overview

jInv consists of five submodules:

1. `ForwardShare` - methods for solving forward problems in parallel.
2. `InverseSolve` - methods commonly used in inverse problems such as misfit functions, regularization and numerical optimization. 
3. `Mesh` - regular and tensor meshes in 2D and 3D as well as differential operators.
4. `LinearSolvers` - interfaces to sparse and (if installed) direct linear solvers that can be used for solving the discretized PDEs.
5. `Utils` - utility functions

# Requirements

jInv is inteded for use with Julia versions 0.4.x.

1. [`KrylovMethods.jl`](https://github.com/lruthotto/KrylovMethods.jl)  - iterative methods for solving (sparse) linear systems. 

Additional (optional) packages for higher performance. `jInv` detects automatically if these packages are installed and uses them by default.

1. [`MUMPS.jl`](https://github.com/JuliaSparse/MUMPS.jl) - wrapper for MUMPS. Used as a direct PDE solver. 
2. [`ParSpMatVec.jl`](https://github.com/lruthotto/ParSpMatVec.jl) - shared memory implementation for sparse matrix vector products.


# Installation

In julia type:
```
Pkg.clone("https://github.com/JuliaInv/jInv.jl","jInv")
Pkg.test("jInv")
```

# Packages using jInv

1. [`DivSigGrad.jl`](https://github.com/JuliaInv/DivSigGrad.jl) - Inverse conductivity problems
2. [`FWI.jl`](https://github.com/JuliaInv/FWI.jl) - Full Waveform Inversion
