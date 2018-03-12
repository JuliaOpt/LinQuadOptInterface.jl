MOI.copy!(dest::LinQuadOptimizer, src::MOI.ModelLike) = MOIU.defaultcopy!(dest, src)
MOI.copy!(dest::MOI.ModelLike, src::LinQuadOptimizer) = MOIU.defaultcopy!(dest, src)
MOI.copy!(dest::LinQuadOptimizer, src::LinQuadOptimizer) = MOIU.defaultcopy!(dest, src)

MOI.get(::LinQuadOptimizer, ::MOI.ListOfVariableAttributesSet) = MOI.AbstractVariableAttribute[]
MOI.get(::LinQuadOptimizer, ::MOI.ListOfModelAttributesSet) = [MOI.ObjectiveSense, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}]
MOI.get(::LinQuadOptimizer, ::MOI.ListOfConstraintAttributesSet{F,S}) where F where S = MOI.AbstractConstraintAttribute[]

