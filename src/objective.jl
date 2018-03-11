function MOI.canset(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{F}) where F<:MOI.AbstractFunction
    return F in lqs_supported_objectives(m)
end

#=
    Set the objective
=#
function MOI.set!(m::LinQuadOptimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    _setsense!(m, sense)
    nothing
end
function MOI.set!(m::LinQuadOptimizer, ::MOI.ObjectiveFunction, objf::Linear)
    cobjf = MOIU.canonical(objf)
    unsafe_set!(m, MOI.ObjectiveFunction{Linear}(), cobjf)
end
function unsafe_set!(m::LinQuadOptimizer, ::MOI.ObjectiveFunction, objf::Linear)
    if m.obj_is_quad
        # previous objective was quadratic...
        m.obj_is_quad = false
        # zero quadratic part
        lqs_copyquad!(m, Int[], Int[], Float64[])
    end
    lqs_chgobj!(m, getcol.(m, objf.variables), objf.coefficients)
    m.objective_constant = objf.constant
    nothing
end

#=
    Set the objective sense
=#

function _setsense!(m::LinQuadOptimizer, sense::MOI.OptimizationSense)
    if sense == MOI.MinSense
        lqs_chgobjsen!(m, :Min)
        m.obj_sense = MOI.MinSense
    elseif sense == MOI.MaxSense
        lqs_chgobjsen!(m, :Max)
        m.obj_sense = MOI.MaxSense
    elseif sense == MOI.FeasibilitySense
        lqs_chgobjsen!(m, :Min)
        unsafe_set!(m, MOI.ObjectiveFunction{Linear}(), MOI.ScalarAffineFunction(VarInd[],Float64[],0.0))
        m.obj_is_quad = false
        m.obj_sense = MOI.FeasibilitySense
    else
        error("Sense $(sense) unknown.")
    end
end

#=
    Get the objective sense
=#

MOI.get(m::LinQuadOptimizer,::MOI.ObjectiveSense) = m.obj_sense#lqs_getobjsen(m)
MOI.canget(m::LinQuadOptimizer, ::MOI.ObjectiveSense) = true

#=
    Get the objective function
=#

function MOI.get(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{Linear})
    variable_coefficients = lqs_getobj(m)
    Linear(m.variable_references, variable_coefficients, m.objective_constant)
end
# can't get quadratic objective functions
MOI.canget(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{S}) where S = false
MOI.canget(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{Linear}) = !m.obj_is_quad

#=
    Modify objective function
=#

function MOI.modifyobjective!(m::LinQuadOptimizer, chg::MOI.ScalarCoefficientChange{Float64})
    col = m.variable_mapping[chg.variable]
    # 0 row is the objective
    lqs_chgcoef!(m, 0, col, chg.new_coefficient)
end
MOI.canmodifyobjective(m::LinQuadOptimizer, ::Type{MOI.ScalarCoefficientChange{Float64}}) = true

#=
    Set quadratic objective
=#

function MOI.set!(m::LinQuadOptimizer, ::MOI.ObjectiveFunction, objf::Quad)
    m.obj_is_quad = true
    lqs_chgobj!(m,
        getcol.(m, objf.affine_variables),
        objf.affine_coefficients
    )
    ri, ci, vi = reduceduplicates(
        getcol.(m, objf.quadratic_rowvariables),
        getcol.(m, objf.quadratic_colvariables),
        objf.quadratic_coefficients
    )
    lqs_copyquad!(m,
        ri,
        ci,
        vi
    )
    m.objective_constant = objf.constant
    nothing
end