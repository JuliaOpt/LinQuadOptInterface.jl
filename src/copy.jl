MOI.copy!(dest::LinQuadOptimizer, src::MOI.ModelLike) = MOIU.defaultcopy!(dest, src)
MOI.copy!(dest::MOI.ModelLike, src::LinQuadOptimizer) = MOIU.defaultcopy!(dest, src)
MOI.copy!(dest::LinQuadOptimizer, src::LinQuadOptimizer) = MOIU.defaultcopy!(dest, src)

MOI.canget(::LinQuadOptimizer, ::MOI.ListOfVariableAttributesSet) = true
MOI.get(::LinQuadOptimizer, ::MOI.ListOfVariableAttributesSet) = MOI.AbstractVariableAttribute[]
MOI.canget(::LinQuadOptimizer, ::MOI.ListOfModelAttributesSet) = true
MOI.get(::LinQuadOptimizer, ::MOI.ListOfModelAttributesSet) = [MOI.ObjectiveSense, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}]
MOI.canget(::LinQuadOptimizer, ::MOI.ListOfConstraintAttributesSet{F,S}) where F where S = true
MOI.get(::LinQuadOptimizer, ::MOI.ListOfConstraintAttributesSet{F,S}) where F where S = MOI.AbstractConstraintAttribute[]

