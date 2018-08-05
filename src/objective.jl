#=
    The Objective Sense
=#

MOI.supports(::LinQuadOptimizer, ::MOI.ObjectiveSense) = true
MOI.canget(::LinQuadOptimizer, ::MOI.ObjectiveSense) = true
MOI.get(model::LinQuadOptimizer,::MOI.ObjectiveSense) = model.obj_sense
function MOI.set!(model::LinQuadOptimizer, ::MOI.ObjectiveSense,
                  sense::MOI.OptimizationSense)
    if sense == MOI.MinSense
        change_objective_sense!(model, :min)
        model.obj_sense = MOI.MinSense
    elseif sense == MOI.MaxSense
        change_objective_sense!(model, :max)
        model.obj_sense = MOI.MaxSense
    elseif sense == MOI.FeasibilitySense
        # we set the objective sense to :min, and the objective to 0.0
        change_objective_sense!(model, :min)
        unsafe_set!(model, MOI.ObjectiveFunction{Linear}(),
                    MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Float64}[],
                    0.0))
        model.obj_type = AffineObjective
        model.obj_sense = MOI.FeasibilitySense
    else
        throw(MOI.CannotSetAttribute(MOI.ObjectiveSense,
                                     "ObjectiveSense $(sense) not recognised."))
    end
end

#=
    The Objective Function
=#

function __assert_objective__(model::LinQuadOptimizer,
                              attribute::MOI.ObjectiveFunction{F}) where F
    if MOI.get(model, MOI.ObjectiveSense()) == MOI.FeasibilitySense
        # it doesn't make sense to set an objective for a feasibility problem
        throw(MOI.CannotSetAttribute(attribute, "Cannot set $(attribute) when" *
            " MOI.ObjectiveSense is MOI.FeasibilitySense."))
    elseif !(F in supported_objectives(model))
        throw(MOI.UnsupportedAttribute(attribute))
    end
end

function MOI.set!(model::LinQuadOptimizer,
                  attribute::MOI.ObjectiveFunction{MOI.SingleVariable},
                  objective::MOI.SingleVariable)
     __assert_objective__(model, attribute)
    if model.obj_type == QuadraticObjective
        set_quadratic_objective!(model, Int[], Int[], Float64[])
    end
    model.obj_type = SingleVariableObjective
    model.single_obj_var = objective.variable
    set_linear_objective!(model, [get_column(model, objective.variable)], [1.0])
    set_constant_objective!(model, 0.0)
end

function MOI.set!(model::LinQuadOptimizer, attribute::MOI.ObjectiveFunction{F},
                  objective::Linear) where F
    __assert_objective__(model, attribute)
    unsafe_set!(model, MOI.ObjectiveFunction{Linear}(), MOIU.canonical(objective))
end

"""
    unsafe_set!(m, ::MOI.ObjectiveFunction{F}, objective::Linear) where F

Sets a linear objective function without cannonicalizing `objective`.
"""
function unsafe_set!(model::LinQuadOptimizer, ::MOI.ObjectiveFunction{F},
                     objective::Linear) where F
    if model.obj_type == QuadraticObjective
        # previous objective was quadratic, so zero quadratic part
        set_quadratic_objective!(model, Int[], Int[], Float64[])
    end
    model.obj_type = AffineObjective
    model.single_obj_var = nothing
    set_linear_objective!(model,
        map(term -> get_column(model, term.variable_index), objective.terms),
        map(term -> term.coefficient, objective.terms)
    )
    set_constant_objective!(model, objective.constant)
end

function MOI.set!(model::LinQuadOptimizer, attribute::MOI.ObjectiveFunction,
                  objective::Quad)
    __assert_objective__(model, attribute)
    model.obj_type = QuadraticObjective
    model.single_obj_var = nothing
    set_linear_objective!(model,
        map(term -> get_column(model, term.variable_index), objective.affine_terms),
        map(term -> term.coefficient, objective.affine_terms)
    )
    columns_1, columns_2, coefficients = reduce_duplicates!(
        map(term -> get_column(model, term.variable_index_1), objective.quadratic_terms),
        map(term -> get_column(model, term.variable_index_2), objective.quadratic_terms),
        map(term -> term.coefficient, objective.quadratic_terms)
    )
    set_quadratic_objective!(model, columns_1, columns_2, coefficients)
    set_constant_objective!(model, objective.constant)
end

#=
    Get the objective function
=#
function MOI.supports(model::LinQuadOptimizer, ::MOI.ObjectiveFunction{F}) where F
    return F in supported_objectives(model)
end

function MOI.canget(model::LinQuadOptimizer, ::MOI.ObjectiveFunction{MOI.SingleVariable})
    return model.obj_type == SingleVariableObjective
end

function MOI.get(model::LinQuadOptimizer, ::MOI.ObjectiveFunction{MOI.SingleVariable})
    if model.obj_type != MOI.SingleVariable
        error("Cannot get SingleVariable objective function as the model has " *
              " a $(model.obj_type) objective function.")
    end
    return MOI.SingleVariable(model.single_obj_var::MOI.VariableIndex)
end

function MOI.canget(model::LinQuadOptimizer, ::MOI.ObjectiveFunction{Linear})
    return model.obj_type != QuadraticObjective
end

function MOI.get(model::LinQuadOptimizer, ::MOI.ObjectiveFunction{Linear})
    if model.obj_type == QuadraticObjective
        error("Cannot get ScalarAffine objective function as the model has " *
              " a quadratic objective function.")
    end
    variable_coefficients = zeros(length(model.variable_references))
    get_linear_objective!(model, variable_coefficients)
    terms = map(
        (variable, coefficient) -> MOI.ScalarAffineTerm{Float64}(coefficient, variable),
        model.variable_references,
        variable_coefficients
    )
    return Linear(terms, get_constant_objective(model))
end

function MOI.canget(model::LinQuadOptimizer, ::MOI.ObjectiveFunction{Quad})
    return model.obj_type == QuadraticObjective
end

function MOI.get(model::LinQuadOptimizer, ::MOI.ObjectiveFunction{Quad})
    if model.obj_type != QuadraticObjective
        error("Cannot get ScalarQuadratic objective function as the model has " *
              " a $(model.obj_type) objective function.")
    end
    variable_coefficients = zeros(length(model.variable_references))
    get_linear_objective!(model, variable_coefficients)
    affine_terms = map(
        (variable, coefficient) -> MOI.ScalarAffineTerm{Float64}(coefficient, variable),
        model.variable_references,
        variable_coefficients
    )
    Q = get_quadratic_terms_objective(model)
    rows = rowvals(Q)
    coefficients = nonzeros(Q)
    quadratic_terms = MOI.ScalarQuadraticTerm{Float64}[]
    sizehint!(quadratic_terms, length(coefficients))
    for (column, variable) in enumerate(model.variable_references)
        for j in nzrange(Q, column)
            row = rows[j]
            push!(quadratic_terms,
                  MOI.ScalarQuadraticTerm{Float64}(
                      coefficients[j],
                      model.variable_references[row],
                      variable)
            )
        end
    end
    return Quad(affine_terms, quadratic_terms, get_constant_objective(model))
end

#=
    Modify objective function
=#

function MOI.modify!(model::LinQuadOptimizer, ::MOI.ObjectiveFunction{F},
                     change::MOI.ScalarCoefficientChange{Float64}) where F<:MOI.AbstractScalarFunction
    if F <: MOI.ScalarQuadraticFunction && model.obj_type != QuadraticObjective
        throw(MOI.UnsupportedObjectiveModification(change,
            "ObjectiveFunction is not a ScalarQuadraticFunction."))
    elseif F <: MOI.ScalarAffineFunction && model.obj_type != AffineObjective
        throw(MOI.UnsupportedObjectiveModification(change,
            "ObjectiveFunction is not a ScalarAffineFunction."))
    end
    if model.obj_type == SingleVariableObjective
        model.obj_type = AffineObjective
        model.single_obj_var = nothing
    end
    change_objective_coefficient!(model, get_column(model, change.variable),
                                  change.new_coefficient)
end

function MOI.modify!(model::LinQuadOptimizer, ::MOI.ObjectiveFunction{F},
                     change::MOI.ScalarConstantChange{Float64}) where F<:MOI.AbstractScalarFunction
    if F == MOI.SingleVariable
        throw(MOI.UnsupportedObjectiveModification(change,
            "ObjectiveFunction is a SingleVariable. Cannot change constant term."))
    end
    set_constant_objective!(model, change.new_constant)
end
