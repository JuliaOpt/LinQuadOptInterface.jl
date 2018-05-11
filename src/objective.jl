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
        unsafe_set!(m, MOI.ObjectiveFunction{Linear}(), MOI.ScalarAffineFunction(VarInd[],Float64[],0.0))
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
    m.objective_constant = 0
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
    set_linear_objective!(m, getcol.(m, objf.variables), objf.coefficients)
    m.objective_constant = objf.constant
    nothing
end

function MOI.set!(m::LinQuadOptimizer, ::MOI.ObjectiveFunction, objf::Quad)
    m.obj_type = QuadraticObjective
    m.single_obj_var = nothing
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

MOI.canget(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{MOI.SingleVariable}) = m.obj_type == SingleVariableObjective
function MOI.get(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{MOI.SingleVariable})
    SingleVariable(get(m.single_obj_var))
end

MOI.canget(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{Linear}) = m.obj_type != QuadraticObjective
function MOI.get(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{Linear})
    variable_coefficients = zeros(length(m.variable_references))
    get_linear_objective!(m, variable_coefficients)
    Linear(m.variable_references, variable_coefficients, m.objective_constant)
end

#=
    Modify objective function
=#

MOI.canmodifyobjective(m::LinQuadOptimizer, ::Type{MOI.ScalarCoefficientChange{Float64}}) = true
function MOI.modifyobjective!(m::LinQuadOptimizer, chg::MOI.ScalarCoefficientChange{Float64})
    if m.obj_type == SingleVariableObjective
        m.obj_type = AffineObjective
        m.single_obj_var = nothing
    end
    col = m.variable_mapping[chg.variable]
    # 0 row is the objective
    change_coefficient!(m, 0, col, chg.new_coefficient)
end
