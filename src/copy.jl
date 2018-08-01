MOI.copy!(dest::LinQuadOptimizer, src::MOI.ModelLike; copynames=false)    = MOIU.defaultcopy!(dest, src, copynames)
MOI.copy!(dest::MOI.ModelLike,    src::LinQuadOptimizer; copynames=false) = MOIU.defaultcopy!(dest, src, copynames)
MOI.copy!(dest::LinQuadOptimizer, src::LinQuadOptimizer; copynames=false) = MOIU.defaultcopy!(dest, src, copynames)

MOI.canget(::LinQuadOptimizer, ::MOI.ListOfVariableAttributesSet) = true
MOI.get(::LinQuadOptimizer, ::MOI.ListOfVariableAttributesSet) = MOI.AbstractVariableAttribute[]

MOI.canget(::LinQuadOptimizer, ::MOI.ListOfModelAttributesSet) = true
MOI.get(m::LinQuadOptimizer, ::MOI.ListOfModelAttributesSet) = [
    MOI.ObjectiveSense(),
    MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
    (MOI.get(m, MOI.Name()) != "") ? MOI.Name() : nothing
]

MOI.canget(::LinQuadOptimizer, ::MOI.ListOfConstraintAttributesSet{F,S}) where F where S = true
MOI.get(::LinQuadOptimizer, ::MOI.ListOfConstraintAttributesSet{F,S}) where F where S = MOI.AbstractConstraintAttribute[]
