## Setup
import Plots
Plots.pyplot();

N = 100;
nSteps = 1000;
## Native
include("./Native/lidCavity.jl");
Main.lidCavity.main(N,100., 1e-5, 1, false);
@time Main.lidCavity.main(N,100., 1e-5, nSteps, false);

## Built-in
include("./Builtin/lidCavity.jl");
Main.lidCavity.main(N,100., 1e-5, 1, false);
@time Main.lidCavity.main(N,100., 1e-5, nSteps, false);

## SIMD
include("./SIMD/lidCavity.jl");
Main.lidCavity.main(N,100., 1e-5, 1, false);
@time Main.lidCavity.main(N,100., 1e-5, nSteps, false);

## Threads
include("./Threads/lidCavity.jl");
Main.lidCavity.main(N,100., 1e-5, 1, false);
@time Main.lidCavity.main(N,100., 1e-5, nSteps, false);

## Broadcasting
#include("./Vectorization/lidCavity.jl");
# Main.lidCavity.main(N,100., 1e-5, 1, false);
#@time Main.lidCavity.main(N,100., 1e-5, nSteps, false);

## Linear Algebra
include("./LinearAlgebra/lidCavity.jl");
Main.lidCavity.main(N,100., "CG", 1e-5, 1, false);
@time Main.lidCavity.main(N,100.,"CG", 1e-5, nSteps, false);

import IterativeSolvers
import SparseArrays
