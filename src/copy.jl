MOI.copy!(dest::LinQuadOptimizer, src::MOI.ModelLike) = MOIU.defaultcopy!(dest, src)
MOI.copy!(dest::MOI.ModelLike, src::LinQuadOptimizer) = MOIU.defaultcopy!(dest, src)
MOI.copy!(dest::LinQuadOptimizer, src::LinQuadOptimizer) = MOIU.defaultcopy!(dest, src)