# jInv - A Flexible Flexible Julia Package for PDE Parameter Estimation

```@contents
	Pages =["Mesh.md","LinearSolvers.md","InverseSolve.md","ForwardShare.md","Utils.md"]
	Depth=2
```

`jInv` is a flexible framework for PDE parameter estimation in Julia. It provides easy to extend core functions used in PDE-constrained inverse problems.

## Build status
[![Build Status](https://travis-ci.org/JuliaInv/jInv.jl.svg?branch=master)](https://travis-ci.org/JuliaInv/jInv.jl) 
[![Coverage Status](https://coveralls.io/repos/github/JuliaInv/jInv.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaInv/jInv.jl?branch=master)
[![Build status](https://ci.appveyor.com/api/projects/status/0pxgtmm08b0w6wgh?svg=true)](https://ci.appveyor.com/project/JuliaInv/jinv-jl-81lel)

# Overview

jInv consists of five submodules:

1. `ForwardShare` - methods for solving forward problems in parallel.
2. `InverseSolve` - methods commonly used in inverse problems such as misfit functions, regularization and numerical optimization. 
3. `Mesh` - regular and tensor meshes in 2D and 3D as well as differential operators.
4. `LinearSolvers` - interfaces to sparse and (if installed) direct linear solvers that can be used for solving the discretized PDEs.
5. `Utils` - utility functions

# Installation

In julia type:

```@repl
Pkg.clone("https://github.com/JuliaInv/jInv.jl","jInv")
Pkg.test("jInv")
```

# Requirements

jInv is developed using Julia versions 0.4.x.

1. [`KrylovMethods.jl`](https://github.com/lruthotto/KrylovMethods.jl)  - iterative methods for solving (sparse) linear systems. 

Additional (optional) packages for higher performance. `jInv` detects automatically if these packages are installed and uses them by default.

1. [`MUMPS.jl`](https://github.com/JuliaSparse/MUMPS.jl) - wrapper for MUMPS. Used as a direct PDE solver. 
2. [`ParSpMatVec.jl`](https://github.com/lruthotto/ParSpMatVec.jl) - shared memory implementation for sparse matrix vector products.


