#=
    The Objective Sense
=#

MOI.canget(m::LinQuadOptimizer, ::MOI.ObjectiveSense) = true
MOI.get(m::LinQuadOptimizer,::MOI.ObjectiveSense) = m.obj_sense

MOI.canset(::LinQuadOptimizer, ::MOI.ObjectiveSense) = true
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

function MOI.canset(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{F}) where F<:MOI.AbstractFunction
    if MOI.get(m, MOI.ObjectiveSense()) == MOI.FeasibilitySense
        # it doesn't make sense to set an objective for a feasibility problem
        return false
    end
    return F in supported_objectives(m)
end

function MOI.set!(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{MOI.SingleVariable}, var::MOI.SingleVariable)
    if m.obj_type == QuadraticObjective
        set_quadratic_objective!(m, Int[], Int[], Float64[])
    end
    m.obj_type = SingleVariableObjective
    m.single_obj_var = var.variable
    set_linear_objective!(m, [getcol(m, var.variable)], [1.0])
    set_constant_objective!(m, 0.0)
end

function MOI.set!(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{F}, objf::Linear) where F
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

function MOI.set!(m::LinQuadOptimizer, ::MOI.ObjectiveFunction, objf::Quad)
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

MOI.canget(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{MOI.SingleVariable}) = m.obj_type == SingleVariableObjective
function MOI.get(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{MOI.SingleVariable})
    SingleVariable(m.single_obj_var::MOI.VariableIndex)
end

MOI.canget(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{Linear}) = m.obj_type != QuadraticObjective
function MOI.get(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{Linear})
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

function MOI.canmodify(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{F}, ::Type{MOI.ScalarCoefficientChange{Float64}}) where F
    if F <: MOI.ScalarQuadraticFunction
        return m.obj_type == QuadraticObjective
    elseif F <: MOI.ScalarAffineFunction
        return m.obj_type == AffineObjective
    end
    return false
end
function MOI.modify!(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{F}, chg::MOI.ScalarCoefficientChange{Float64}) where F
    if m.obj_type == SingleVariableObjective
        m.obj_type = AffineObjective
        m.single_obj_var = nothing
    end
    col = m.variable_mapping[chg.variable]
    change_objective_coefficient!(m, col, chg.new_coefficient)
end

function MOI.canmodify(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{F}, ::Type{MOI.ScalarConstantChange{Float64}}) where F
    return !(F == MOI.SingleVariable)
end
function MOI.modify!(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{F}, chg::MOI.ScalarConstantChange{Float64}) where F
    set_constant_objective!(m, chg.new_constant)
end
