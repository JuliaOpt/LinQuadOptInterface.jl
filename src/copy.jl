function MOI.copy!(dest::LinQuadOptimizer, src; copynames=false)
    src::MOI.ModelLike
    return MOIU.defaultcopy!(dest, src, copynames)
end

function MOI.copy!(dest::LinQuadOptimizer, src::LinQuadOptimizer; copynames=false)
    return MOIU.defaultcopy!(dest, src, copynames)
end

MOI.canget(::LinQuadOptimizer, ::MOI.ListOfVariableAttributesSet) = true
function MOI.get(::LinQuadOptimizer, ::MOI.ListOfVariableAttributesSet)
    return MOI.AbstractVariableAttribute[]
end

MOI.canget(::LinQuadOptimizer, ::MOI.ListOfModelAttributesSet) = true
function MOI.get(model::LinQuadOptimizer, ::MOI.ListOfModelAttributesSet)
    return [MOI.ObjectiveSense(),
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        # TODO(odow): why is `nothing` returned?
        (MOI.get(model, MOI.Name()) != "") ? MOI.Name() : nothing
    ]
end

function MOI.canget(::LinQuadOptimizer, ::MOI.ListOfConstraintAttributesSet{F,S}) where {F, S}
    return true
end
function MOI.get(::LinQuadOptimizer, ::MOI.ListOfConstraintAttributesSet{F,S}) where {F, S}
    return MOI.AbstractConstraintAttribute[]
end
