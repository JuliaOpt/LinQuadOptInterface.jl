#=
    Set the objective
=#

function MOI.setobjective!(m::LinQuadSolverInstance, sense::MOI.OptimizationSense, objf::Linear)
    if m.obj_is_quad
        # previous objective was quadratic...
        m.obj_is_quad = false
        # zero quadratic part
        lqs_copyquad!(m.inner, Int[], Int[], Float64[])
    end
    lqs_chgobj!(m.inner, getcol.(m, objf.variables), objf.coefficients)
    _setsense!(m, sense)
    m.objective_constant = objf.constant
    nothing
end

#=
    Set the objective sense
=#

function _setsense!(m::LinQuadSolverInstance, sense::MOI.OptimizationSense)
    if sense == MOI.MinSense
        lqs_chgobjsen!(m.inner, :Min)
    elseif sense == MOI.MaxSense
        lqs_chgobjsen!(m.inner, :Max)
    elseif sense == MOI.FeasibilitySense
        warn("FeasibilitySense not supported. Using MinSense")
        lqs_chgobjsen!(m.inner, :Min)
    else
        error("Sense $(sense) unknown.")
    end
end

#=
    Get the objective sense
=#

MOI.getattribute(m::LinQuadSolverInstance,::MOI.ObjectiveSense) = lqs_getobjsen(m.inner)
MOI.cangetattribute(m::LinQuadSolverInstance, ::MOI.ObjectiveSense) = true

#=
    Get the objective function
=#

function MOI.getattribute(m::LinQuadSolverInstance, ::MOI.ObjectiveFunction)
    variable_coefficients = lqs_getobj(m.inner)
    MOI.ScalarAffineFunction(m.variable_references, variable_coefficients, m.objective_constant)
end
# can't get quadratic objective functions
MOI.cangetattribute(m::LinQuadSolverInstance, ::MOI.ObjectiveFunction) = !m.obj_is_quad

#=
    Modify objective function
=#

function MOI.modifyobjective!(m::LinQuadSolverInstance, chg::MOI.ScalarCoefficientChange{Float64})
    col = m.variable_mapping[chg.variable]
    # 0 row is the objective
    lqs_chgcoef!(m.inner, 0, col, chg.new_coefficient)
end
MOI.canmodifyobjective(m::LinQuadSolverInstance, chg::MOI.ScalarCoefficientChange{Float64}) = true

#=
    Set quadratic objective
=#

function MOI.setobjective!(m::LinQuadSolverInstance, sense::MOI.OptimizationSense, objf::Quad)
    m.obj_is_quad = true
    lqs_chgobj!(m.inner,
        getcol.(m, objf.affine_variables),
        objf.affine_coefficients
    )
    ri, ci, vi = reduceduplicates(
        getcol.(m, objf.quadratic_rowvariables),
        getcol.(m, objf.quadratic_colvariables),
        objf.quadratic_coefficients
    )
    lqs_copyquad!(m.inner,
        ri,
        ci,
        vi
    )
    _setsense!(m, sense)
    m.objective_constant = objf.constant
    nothing
end