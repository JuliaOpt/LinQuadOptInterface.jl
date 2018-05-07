#=
    The Objective Sense
=#

MOI.canget(m::LinQuadOptimizer, ::MOI.ObjectiveSense) = true
MOI.get(m::LinQuadOptimizer,::MOI.ObjectiveSense) = m.obj_sense

MOI.canset(::LinQuadOptimizer, ::MOI.ObjectiveSense) = true
function MOI.set!(m::LinQuadOptimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    if sense == MOI.MinSense
        change_objectivesense!(m, :min)
        m.obj_sense = MOI.MinSense
    elseif sense == MOI.MaxSense
        change_objectivesense!(m, :max)
        m.obj_sense = MOI.MaxSense
    elseif sense == MOI.FeasibilitySense
        # we set the objective sense to :min, and the objective to 0.0
        change_objectivesense!(m, :min)
        unsafe_set!(m, MOI.ObjectiveFunction{Linear}(), MOI.ScalarAffineFunction(VarInd[],Float64[],0.0))
        m.obj_is_quad = false
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
    return F in lqs_supported_objectives(m)
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
    if m.obj_is_quad
        # previous objective was quadratic...
        m.obj_is_quad = false
        # zero quadratic part
        set_quadratic_objective!(m, Int[], Int[], Float64[])
    end
    set_linear_objective!(m, getcol.(m, objf.variables), objf.coefficients)
    m.objective_constant = objf.constant
    nothing
end

function MOI.set!(m::LinQuadOptimizer, ::MOI.ObjectiveFunction, objf::Quad)
    m.obj_is_quad = true
    set_linear_objective!(m,
        getcol.(m, objf.affine_variables),
        objf.affine_coefficients
    )
    ri, ci, vi = reduceduplicates(
        getcol.(m, objf.quadratic_rowvariables),
        getcol.(m, objf.quadratic_colvariables),
        objf.quadratic_coefficients
    )
    set_quadratic_objective!(m,
        ri,
        ci,
        vi
    )
    m.objective_constant = objf.constant
    nothing
end

#=
    Get the objective function
=#

MOI.canget(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{Linear}) = !m.obj_is_quad
function MOI.get(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{Linear})
    variable_coefficients = get_linearobjective(m)
    Linear(m.variable_references, variable_coefficients, m.objective_constant)
end

#=
    Modify objective function
=#

MOI.canmodifyobjective(m::LinQuadOptimizer, ::Type{MOI.ScalarCoefficientChange{Float64}}) = true
function MOI.modifyobjective!(m::LinQuadOptimizer, chg::MOI.ScalarCoefficientChange{Float64})
    col = m.variable_mapping[chg.variable]
    # 0 row is the objective
    change_coefficient!(m, 0, col, chg.new_coefficient)
end
