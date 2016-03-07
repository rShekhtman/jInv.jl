[![Build Status](https://travis-ci.org/JuliaInv/jInv.jl.svg?branch=master)](https://travis-ci.org/JuliaInv/jInv.jl) 
[![Coverage Status](https://coveralls.io/repos/github/JuliaInv/jInv.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaInv/jInv.jl?branch=master)

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
