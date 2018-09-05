
# this function is only defined with `LinQuadOptimizer` in dest because a version
# with `LinQuadOptimizer` in src and anything in dest does not work due to ambiguity.
function MOI.copy_to(dest::LinQuadOptimizer, src; copynames=false)
    src::MOI.ModelLike
    return MOIU.default_copy_to(dest, src, copynames)
end

function MOI.copy_to(dest::LinQuadOptimizer, src::LinQuadOptimizer; copy_names=false)
    return MOIU.default_copy_to(dest, src, copy_names)
end

function MOI.get(::LinQuadOptimizer, ::MOI.ListOfVariableAttributesSet)
    return MOI.AbstractVariableAttribute[]
end

function MOI.get(model::LinQuadOptimizer, ::MOI.ListOfModelAttributesSet)
    return [MOI.ObjectiveSense(),
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        # TODO(odow): why is `nothing` returned?
        (MOI.get(model, MOI.Name()) != "") ? MOI.Name() : nothing
    ]
end

function MOI.get(::LinQuadOptimizer, ::MOI.ListOfConstraintAttributesSet{F,S}) where {F, S}
    return MOI.AbstractConstraintAttribute[]
end
