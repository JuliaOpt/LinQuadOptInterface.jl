#=
    Set the objective
=#
MOI.canset(::LinQuadOptimizer, ::MOI.ObjectiveSense) = true
MOI.canset(::LinQuadOptimizer, ::Type{MOI.ObjectiveSense}) = true
MOI.set!(m::LinQuadOptimizer, ::Type{MOI.ObjectiveSense}, sense::MOI.OptimizationSense) = MOI.set!(m, MOI.ObjectiveSense(), sense)
function MOI.set!(m::LinQuadOptimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    _setsense!(m, sense)
    nothing
end

MOI.canset(m::LinQuadOptimizer, ::Type{MOI.ObjectiveFunction{F}}) where F<:MOI.AbstractFunction = MOI.canset(m, MOI.ObjectiveFunction{F}())
function MOI.canset(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{F}) where F<:MOI.AbstractFunction
    return F in lqs_supported_objectives(m)
end
MOI.set!(m::LinQuadOptimizer, ::Type{MOI.ObjectiveFunction{F}}, objf::Linear) where F = MOI.set!(m, MOI.ObjectiveFunction{F}(), objf::Linear)
function MOI.set!(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{F}, objf::Linear) where F
    cobjf = MOIU.canonical(objf)
    unsafe_set!(m, MOI.ObjectiveFunction{Linear}(), cobjf)
end
function unsafe_set!(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{F}, objf::Linear) where F
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
        lqs_chgobjsen!(m, :min)
        m.obj_sense = MOI.MinSense
    elseif sense == MOI.MaxSense
        lqs_chgobjsen!(m, :max)
        m.obj_sense = MOI.MaxSense
    elseif sense == MOI.FeasibilitySense
        lqs_chgobjsen!(m, :min)
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
MOI.get(m::LinQuadOptimizer,::Type{MOI.ObjectiveSense}) = m.obj_sense#lqs_getobjsen(m)
MOI.canget(m::LinQuadOptimizer, ::MOI.ObjectiveSense) = true
MOI.canget(m::LinQuadOptimizer, ::Type{MOI.ObjectiveSense}) = true

#=
    Get the objective function
=#

MOI.get(m::LinQuadOptimizer, ::Type{MOI.ObjectiveFunction{Linear}}) = MOI.get(m, MOI.ObjectiveFunction{Linear}())
function MOI.get(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{Linear})
    variable_coefficients = lqs_getobj(m)
    Linear(m.variable_references, variable_coefficients, m.objective_constant)
end
# can't get quadratic objective functions
MOI.canget(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{S}) where S = false
MOI.canget(m::LinQuadOptimizer, ::Type{MOI.ObjectiveFunction{S}}) where S = false
MOI.canget(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{Linear}) = !m.obj_is_quad
MOI.canget(m::LinQuadOptimizer, ::Type{MOI.ObjectiveFunction{Linear}}) = !m.obj_is_quad

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
