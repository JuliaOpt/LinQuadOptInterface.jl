# LinQuadOptInterface.jl

| **Build Status** | **Social** |
|:-----------------:|:----------:|
| [![Build Status][build-img]][build-url] [![Codecov branch][codecov-img]][codecov-url] | [![Gitter][gitter-img]][gitter-url] [<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/a/af/Discourse_logo.png/799px-Discourse_logo.png" width="64">][discourse-url] |

[build-img]: https://travis-ci.org/JuliaOpt/LinQuadOptInterface.jl.svg?branch=master
[build-url]: https://travis-ci.org/JuliaOpt/LinQuadOptInterface.jl
[codecov-img]: http://codecov.io/github/JuliaOpt/LinQuadOptInterface.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/JuliaOpt/LinQuadOptInterface.jl?branch=master

[gitter-url]: https://gitter.im/JuliaOpt/JuMP-dev?utm_source=share-link&utm_medium=link&utm_campaign=share-link
[gitter-img]: https://badges.gitter.im/JuliaOpt/JuMP-dev.svg
[discourse-url]: https://discourse.julialang.org/c/domain/opt

LinQuadOptInterface.jl (LQOI) is designed to provide an intermediate interface
to [MathOptInterface.jl](https://github.com/JuliaOpt/MathOptInterface.jl)
for some solvers. The target use-cases are low-level wrappers designed to bridge
low-level mixed integer linear and quadratic solvers.

Examples of packages currently using LQOI include [Clp.jl](https://github.com/JuliaOpt/Clp.jl),
[GLPK.jl](https://github.com/JuliaOpt/GLPK.jl), [Gurobi.jl](https://github.com/JuliaOpt/Gurobi.jl), and [Xpress.jl](https://github.com/JuliaOpt/Xpress.jl).

The interface is documented [here](https://github.com/JuliaOpt/LinQuadOptInterface.jl/blob/master/src/solver_interface.jl).

## Note to solver developers

The use of LQOI for MOI wrappers of low-level solvers is entirely optional.
Using LQOI introduces an extra abstraction layer between a solver and MOI. We
recommend that you carefully analyze the solver's low-level API to check if it
is close to what LQOI expects.

If a solver low-level API does not support most of the functions required by LQOI, then following the example of
[MathOptInterfaceSCS.jl](https://github.com/JuliaOpt/MathOptInterfaceSCS.jl) and
[MathOptInterfaceECOS.jl](https://github.com/JuliaOpt/MathOptInterfaceECOS.jl)
might be a better idea.
