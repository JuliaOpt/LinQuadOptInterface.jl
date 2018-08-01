#=
    The Objective Sense
=#
MOI.supports(::LinQuadOptimizer, ::MOI.ObjectiveSense) = true
MOI.canget(m::LinQuadOptimizer, ::MOI.ObjectiveSense) = true
MOI.get(m::LinQuadOptimizer,::MOI.ObjectiveSense) = m.obj_sense
function MOI.set!(m::LinQuadOptimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    if sense == MOI.MinSense
        change_objective_sense!(m, :min)
        m.obj_sense = MOI.MinSense
    elseif sense == MOI.MaxSense
        change_objective_sense!(m, :max)
        m.obj_sense = MOI.MaxSense
    elseif sense == MOI.FeasibilitySense
        # we set the objective sense to :min, and the objective to 0.0
        change_objective_sense!(m, :min)
        unsafe_set!(m, MOI.ObjectiveFunction{Linear}(), MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Float64}[],0.0))
        m.obj_type = AffineObjective
        m.obj_sense = MOI.FeasibilitySense
    else
        error("Sense $(sense) unknown.")
    end
end

#=
    The Objective Function
=#

function _assert_objective(model::LinQuadOptimizer, attribute::MOI.ObjectiveFunction{F}) where F
    if MOI.get(model, MOI.ObjectiveSense()) == MOI.FeasibilitySense
        # it doesn't make sense to set an objective for a feasibility problem
        error("Cannot set objective function when MOI.ObjectiveSense is MOI.FeasibilitySense.")
    elseif !(F in supported_objectives(model))
        throw(MOI.UnsupportedAttribute(attribute))
    end
end

function MOI.set!(m::LinQuadOptimizer, attr::MOI.ObjectiveFunction{MOI.SingleVariable}, var::MOI.SingleVariable)
     _assert_objective(m, attr)
    if m.obj_type == QuadraticObjective
        set_quadratic_objective!(m, Int[], Int[], Float64[])
    end
    m.obj_type = SingleVariableObjective
    m.single_obj_var = var.variable
    set_linear_objective!(m, [getcol(m, var.variable)], [1.0])
    set_constant_objective!(m, 0.0)
end

function MOI.set!(m::LinQuadOptimizer, attr::MOI.ObjectiveFunction{F}, objf::Linear) where F
    _assert_objective(m, attr)
    cobjf = MOIU.canonical(objf)
    unsafe_set!(m, MOI.ObjectiveFunction{Linear}(), cobjf)
end

"""
    unsafe_set!(m, ::MOI.ObjectiveFunction{F}, objective::Linear) where F

Sets a linear objective function without cannonicalizing `objective`.
"""
function unsafe_set!(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{F}, objf::Linear) where F
    if m.obj_type == QuadraticObjective
        # previous objective was quadratic...
        # zero quadratic part
        set_quadratic_objective!(m, Int[], Int[], Float64[])
    end
    m.obj_type = AffineObjective
    m.single_obj_var = nothing
    columns = [getcol(m, term.variable_index) for term in objf.terms]
    coefficients = [term.coefficient for term in objf.terms]
    set_linear_objective!(m, columns, coefficients)
    set_constant_objective!(m, objf.constant)
    nothing
end

function MOI.set!(m::LinQuadOptimizer, attr::MOI.ObjectiveFunction, objf::Quad)
    _assert_objective(m, attr)
    m.obj_type = QuadraticObjective
    m.single_obj_var = nothing
    columns = [getcol(m, term.variable_index) for term in objf.affine_terms]
    coefficients = [term.coefficient for term in objf.affine_terms]
    set_linear_objective!(m, columns, coefficients)

    quadratic_columns_1 = [getcol(m, term.variable_index_1) for term in objf.quadratic_terms]
    quadratic_columns_2 = [getcol(m, term.variable_index_2) for term in objf.quadratic_terms]
    quadratic_coefficients = [term.coefficient for term in objf.quadratic_terms]
    ri, ci, vi = reduce_duplicates!(
        quadratic_columns_1,
        quadratic_columns_2,
        quadratic_coefficients
    )
    set_quadratic_objective!(m, ri, ci, vi)
    set_constant_objective!(m, objf.constant)
    nothing
end

#=
    Get the objective function
=#
function MOI.supports(model::LinQuadOptimizer, ::MOI.ObjectiveFunction{F}) where F
    return F in supported_objectives(model)
end

MOI.canget(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{MOI.SingleVariable}) = m.obj_type == SingleVariableObjective

function MOI.get(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{MOI.SingleVariable})
    if m.obj_type != MOI.SingleVariable
        error("Cannot get objective function.")
    end
    return MOI.SingleVariable(m.single_obj_var::MOI.VariableIndex)
end

MOI.canget(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{Linear}) = m.obj_type != QuadraticObjective

function MOI.get(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{Linear})
    if m.obj_type == QuadraticObjective
        error("Cannot get objective function.")
    end
    variable_coefficients = zeros(length(m.variable_references))
    get_linear_objective!(m, variable_coefficients)
    terms = map(
        (v,c)->MOI.ScalarAffineTerm{Float64}(c,v),
        m.variable_references,
        variable_coefficients
    )
    Linear(terms, get_constant_objective(m))
end

MOI.canget(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{Quad}) = m.obj_type == QuadraticObjective

function MOI.get(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{Quad})
    if m.obj_type != QuadraticObjective
        error("Cannot get objective function.")
    end
    variable_coefficients = zeros(length(m.variable_references))
    get_linear_objective!(m, variable_coefficients)
    affine_terms = map(
        (v,c)->MOI.ScalarAffineTerm{Float64}(c,v),
        m.variable_references,
        variable_coefficients
    )
    Q = get_quadratic_terms_objective(m)
    rows = rowvals(Q)
    vals = nonzeros(Q)
    nrows, ncols = size(Q)
    quadratic_terms = MOI.ScalarQuadraticTerm{Float64}[]
    sizehint!(quadratic_terms, length(vals))
    for i = 1:ncols
        for j in nzrange(Q, i)
            row = rows[j]
            val = vals[j]
            push!(quadratic_terms, MOI.ScalarQuadraticTerm{Float64}(val, m.variable_references[row], m.variable_references[i]))
        end
    end
    Quad(affine_terms, quadratic_terms, get_constant_objective(m))
end

#=
    Modify objective function
=#

function MOI.modify!(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{F}, chg::MOI.ScalarCoefficientChange{Float64}) where F<:MOI.AbstractScalarFunction
    if F <: MOI.ScalarQuadraticFunction && m.obj_type != QuadraticObjective
        # TODO(odow): don't we want a descriptive error?
        throw(MOI.UnsupportedObjectiveModification(chg))
    elseif F <: MOI.ScalarAffineFunction && m.obj_type != AffineObjective
        # TODO(odow): don't we want a descriptive error?
        throw(MOI.UnsupportedObjectiveModification(chg))
    end
    if m.obj_type == SingleVariableObjective
        m.obj_type = AffineObjective
        m.single_obj_var = nothing
    end
    col = m.variable_mapping[chg.variable]
    change_objective_coefficient!(m, col, chg.new_coefficient)
end

function MOI.modify!(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{F}, chg::MOI.ScalarConstantChange{Float64}) where F<:MOI.AbstractScalarFunction
    if F == MOI.SingleVariable
        # TODO(odow): don't we want a descriptive error?
        throw(MOI.UnsupportedObjectiveModification(chg))
    end
    set_constant_objective!(m, chg.new_constant)
end
