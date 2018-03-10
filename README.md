# LinQuadOptInterface.jl

LinQuadOptInterface.jl (LQOI) is designed to simplify [MathOptInterface.jl](https://github.com/JuliaOpt/MathOptInterface.jl)'s (MOI) implementation. The target use cases are low-level wrappers designed to bridge low-level Integer Linear and Quadratic solvers, for instance [GLPK.jl](https://github.com/JuliaOpt/GLPK.jl), [Gurobi.jl](https://github.com/JuliaOpt/Gurobi.jl), [Xpress.jl](https://github.com/JuliaOpt/Xpress.jl) and [CPLEX.jl](https://github.com/JuliaOpt/CPLEX.jl).

The use of LQOI for MOI implementations is entirely optional. Using LQOI introduces an extra abstraction layer between a solver and MOI. Its recommended to carefully analyse if the solver's low-level API is close to what LQOI expects, otherwise a direct implementation of MOI might be a better option.

## <a name="modifications"></a> Problem Modifications

LQOI provides access to many MOI functionalities mainly related to problem modifications:

1. add constraints/variables 1-by-1 and in batches
2. remove constraints/variables
3. change constraint coefficients
4. change constraint type
5. change constraint rhs
6. change variable bounds and type

## LinQuadOptimizer

In LinQuadOptInterface.jl the MOI [`AbstractOptimizer`](http://www.juliaopt.org/MathOptInterface.jl/latest/apireference.html#MathOptInterface.AbstractOptimizer) is specialized to [`LinQuadOptimizer`](https://github.com/JuliaOpt/LinQuadOptInterface.jl/blob/99b2a3dfe78e000330475f08766f6681ecf633ab/src/LinQuadOptInterface.jl#L131). In this implementation all the above mentioned modifications are carried out by the low-level solver's own functionalities and hence the `LinQuadOptimizer` can be used without a [`CachingOptimizer`](https://github.com/JuliaOpt/MathOptInterface.jl/blob/60c5ee85addb65ada33cb1d922691f23e5a518e2/src/Utilities/cachingoptimizer.jl#L8). The latter will keep all problem data in the Julia level and typically push to a low-level solver the complete problem at once.

In contrast, since `LinQuadOptimizer` incrementally pushes data to the low-level solver, it stores a small subset of the problem data at the Julia level; typically only references to constraints and variables.

## Current uses

This package is currently used to implement [MathOptInterfaceXpress.jl](https://github.com/JuliaOpt/MathOptInterfaceXpress.jl), [MathOptInterfaceCPLEX.jl](https://github.com/JuliaOpt/MathOptInterfaceCPLEX.jl), [MathOptInterfaceGurobi.jl](https://github.com/JuliaOpt/MathOptInterfaceGurobi.jl) and [MathOptInterfaceGLPK.jl](https://github.com/JuliaOpt/MathOptInterfaceGLPK.jl). The last one being only a integer linear solver.

All these solvers have low-level APIs which supports most of the above metioned [modifications](#modifications). Hence, data storage is simplified and duplications are avoided.

## Other possible uses

This package is only recommended if a solver low-level API which supports most of the above mentioned [modifications](#modifications).

If a solver low-level API does not support most of the above mentioned [modifications](#modifications), then following the example of [MathOptInterfaceSCS.jl](https://github.com/JuliaOpt/MathOptInterfaceSCS.jl) and [MathOptInterfaceECOS.jl](https://github.com/JuliaOpt/MathOptInterfaceECOS.jl) might be a better idea.
