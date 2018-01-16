#=
    Set the objective
=#
function MOI.set!(m::LinQuadSolverInstance, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    _setsense!(m, sense)
    nothing
end
function MOI.set!(m::LinQuadSolverInstance, ::MOI.ObjectiveFunction, objf::Linear)
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

function _setsense!(m::LinQuadSolverInstance, sense::MOI.OptimizationSense)
    if sense == MOI.MinSense
        lqs_chgobjsen!(m, :Min)
    elseif sense == MOI.MaxSense
        lqs_chgobjsen!(m, :Max)
    elseif sense == MOI.FeasibilitySense
        warn("FeasibilitySense not supported. Using MinSense")
        lqs_chgobjsen!(m, :Min)
    else
        error("Sense $(sense) unknown.")
    end
end

#=
    Get the objective sense
=#

MOI.get(m::LinQuadSolverInstance,::MOI.ObjectiveSense) = lqs_getobjsen(m)
MOI.canget(m::LinQuadSolverInstance, ::MOI.ObjectiveSense) = true

#=
    Get the objective function
=#

function MOI.get(m::LinQuadSolverInstance, ::MOI.ObjectiveFunction)
    variable_coefficients = lqs_getobj(m)
    coefs = Float64[]
    vars = VarInd[]
    for i in eachindex(variable_coefficients)
        if abs(variable_coefficients[i]) > 1e-20
            push!(coefs, variable_coefficients[i])
            push!(vars, m.variable_references[i])
        end
    end
    MOI.ScalarAffineFunction(vars, coefs, m.objective_constant)
end
# can't get quadratic objective functions
MOI.canget(m::LinQuadSolverInstance, ::MOI.ObjectiveFunction) = !m.obj_is_quad

#=
    Modify objective function
=#

function MOI.modifyobjective!(m::LinQuadSolverInstance, chg::MOI.ScalarCoefficientChange{Float64})
    col = m.variable_mapping[chg.variable]
    # 0 row is the objective
    lqs_chgcoef!(m, 0, col, chg.new_coefficient)
end
MOI.canmodifyobjective(m::LinQuadSolverInstance, chg::MOI.ScalarCoefficientChange{Float64}) = true

#=
    Set quadratic objective
=#

function MOI.set!(m::LinQuadSolverInstance, ::MOI.ObjectiveFunction, objf::Quad)
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