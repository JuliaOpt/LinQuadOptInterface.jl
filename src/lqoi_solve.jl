function hasquadratic(m::LinQuadSolverInstance)
    m.obj_is_quad || (length(cmap(m).q_less_than) + length(cmap(m).q_greater_than) + length(cmap(m).q_equal_to) > 0)
end

#=
    Optimize the model
=#

function MOI.optimize!(m::LinQuadSolverInstance)
    # reset storage
    fill!(m.variable_primal_solution, NaN)
    fill!(m.variable_dual_solution, NaN)
    fill!(m.constraint_primal_solution, NaN)
    fill!(m.constraint_dual_solution, NaN)
    m.primal_status = MOI.UnknownResultStatus
    m.dual_status   = MOI.UnknownResultStatus
    m.primal_result_count = 0
    m.dual_result_count = 0

    t = time()
    if hasinteger(m)
        lqs_mipopt!(m.inner)
    elseif hasquadratic(m)
        lqs_qpopt!(m.inner)
    else
        lqs_lpopt!(m.inner)
    end
    m.solvetime = time() - t

    # termination_status
    m.termination_status = lqs_terminationstatus(m)
    m.primal_status = lqs_primalstatus(m)
    m.dual_status = lqs_dualstatus(m)

    if m.primal_status in [MOI.FeasiblePoint, MOI.InfeasiblePoint]
        # primal solution exists
        lqs_getx!(m.inner, m.variable_primal_solution)
        lqs_getax!(m.inner, m.constraint_primal_solution)
        m.primal_result_count = 1
        # CPLEX can return infeasible points
    elseif m.primal_status == MOI.InfeasibilityCertificate
        lqs_getray!(m.inner, m.variable_primal_solution)
        m.primal_result_count = 1
    end
    if m.dual_status in [MOI.FeasiblePoint, MOI.InfeasiblePoint]
        # dual solution exists
        lqs_getdj!(m.inner, m.variable_dual_solution)
        lqs_getpi!(m.inner, m.constraint_dual_solution)
        m.dual_result_count = 1
        # dual solution may not be feasible
    elseif m.dual_status == MOI.InfeasibilityCertificate
        lqs_dualfarkas!(m.inner, m.constraint_dual_solution)
        m.dual_result_count = 1
    end

    #=
        CPLEX has the dual convention that the sign of the dual depends on the
        optimization sense. This isn't the same as the MOI convention so we need
        to correct that.
    =#
    # TODO
    if MOI.getattribute(m, MOI.ObjectiveSense()) == MOI.MaxSense
        m.constraint_dual_solution *= -1
        m.variable_dual_solution *= -1
    end
@show pwd()
    MOI.writeproblem(m, "debug.lp", "l")
    @show m.variable_primal_solution
end


#=
    Result Count
=#
function MOI.getattribute(m::LinQuadSolverInstance, ::MOI.ResultCount)
    max(m.primal_result_count, m.dual_result_count)
end
MOI.cangetattribute(m::LinQuadSolverInstance, ::MOI.ResultCount) = true

#=
    Termination status
=#

function MOI.getattribute(m::LinQuadSolverInstance, ::MOI.TerminationStatus)
    m.termination_status
end
MOI.cangetattribute(m::LinQuadSolverInstance, ::MOI.TerminationStatus) = true

#=
    Primal status
=#

function MOI.getattribute(m::LinQuadSolverInstance, p::MOI.PrimalStatus)
    m.primal_status
end
function MOI.cangetattribute(m::LinQuadSolverInstance, p::MOI.PrimalStatus)
    m.primal_result_count >= p.N
end

#=
    Dual status
=#

function MOI.getattribute(m::LinQuadSolverInstance, d::MOI.DualStatus)
    m.dual_status
end
function MOI.cangetattribute(m::LinQuadSolverInstance, d::MOI.DualStatus)
    m.dual_result_count >= d.N
end

#=
    Objective Value
=#


function MOI.getattribute(m::LinQuadSolverInstance, attr::MOI.ObjectiveValue)
    if attr.resultindex == 1
        lqs_getobjval(m.inner) + m.objective_constant
    else
        error("Unable to access multiple objective values")
    end
end
function MOI.cangetattribute(m::LinQuadSolverInstance, attr::MOI.ObjectiveValue)
    if attr.resultindex == 1
        return true
    else
        return false
    end
end

#=
    Variable Primal solution
=#


function MOI.getattribute(m::LinQuadSolverInstance, ::MOI.VariablePrimal, v::MOI.VariableReference)
    col = m.variable_mapping[v]
    return m.variable_primal_solution[col]
end
MOI.cangetattribute(m::LinQuadSolverInstance, ::MOI.VariablePrimal, v::MOI.VariableReference) = true

function MOI.getattribute(m::LinQuadSolverInstance, ::MOI.VariablePrimal, v::Vector{MOI.VariableReference})
    MOI.getattribute.(m, MOI.VariablePrimal(), v)
end
MOI.cangetattribute(m::LinQuadSolverInstance, ::MOI.VariablePrimal, v::Vector{MOI.VariableReference}) = true

#=
    Variable Dual solution
=#


function MOI.getattribute(m::LinQuadSolverInstance,::MOI.ConstraintDual, c::SVCR{<: Union{LE, GE, EQ, IV}})
    vref = m[c]
    col = m.variable_mapping[vref]
    return m.variable_dual_solution[col]
end
MOI.cangetattribute(m::LinQuadSolverInstance, ::MOI.ConstraintDual, c::SVCR{<: Union{LE, GE, EQ, IV}}) = true

#=
    Constraint Primal solution
=#

function MOI.getattribute(m::LinQuadSolverInstance, ::MOI.ConstraintPrimal, c::LCR{<: Union{LE, GE, EQ, IV}})
    row = m[c]
    return m.constraint_primal_solution[row]
end
MOI.cangetattribute(m::LinQuadSolverInstance, ::MOI.ConstraintPrimal,c::LCR{<: Union{LE, GE, EQ, IV}}) = true


# vector valued constraint duals
MOI.getattribute(m::LinQuadSolverInstance, ::MOI.ConstraintPrimal, c::VLCR{<: Union{MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives}}) = m.constraint_primal_solution[m[c]]
MOI.cangetattribute(m::LinQuadSolverInstance, ::MOI.ConstraintPrimal,c::VLCR{<: Union{MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives}}) = true

#=
    Constraint Dual solution
=#

_checkdualsense(::LCR{LE}, dual) = dual <= 0.0
_checkdualsense(::LCR{GE}, dual) = dual >= 0.0
_checkdualsense(::LCR{IV}, dual) = true
_checkdualsense(::LCR{EQ}, dual) = true

function MOI.getattribute(m::LinQuadSolverInstance, ::MOI.ConstraintDual, c::LCR{<: Union{LE, GE, EQ, IV}})
    row = m[c]
    dual = m.constraint_dual_solution[row]
    @assert _checkdualsense(c, dual)
    return dual
end
MOI.cangetattribute(m::LinQuadSolverInstance, ::MOI.ConstraintDual,c::LCR{<: Union{LE, GE, EQ, IV}}) = true

# vector valued constraint duals
MOI.getattribute(m::LinQuadSolverInstance, ::MOI.ConstraintDual, c::VLCR{<: Union{MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives}}) = m.constraint_dual_solution[m[c]]
MOI.cangetattribute(m::LinQuadSolverInstance, ::MOI.ConstraintDual,c::VLCR{<: Union{MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives}}) = true

#=
    Solution Attributes
=#

# struct ObjectiveBound <: AbstractSolverInstanceAttribute end
MOI.getattribute(m::LinQuadSolverInstance, ::MOI.ObjectiveBound) = lqs_getbestobjval(m.inner)
MOI.cangetattribute(m::LinQuadSolverInstance, ::MOI.ObjectiveBound) = true

# struct RelativeGap <: AbstractSolverInstanceAttribute  end
MOI.getattribute(m::LinQuadSolverInstance, ::MOI.RelativeGap) = lqs_getmiprelgap(m.inner)
MOI.cangetattribute(m::LinQuadSolverInstance, ::MOI.RelativeGap) = true

# struct SolveTime <: AbstractSolverInstanceAttribute end
MOI.getattribute(m::LinQuadSolverInstance, ::MOI.SolveTime) = m.solvetime
MOI.cangetattribute(m::LinQuadSolverInstance, ::MOI.SolveTime) = true

# struct SimplexIterations <: AbstractSolverInstanceAttribute end
MOI.getattribute(m::LinQuadSolverInstance, ::MOI.SimplexIterations) = lqs_getitcnt(m.inner)
MOI.cangetattribute(m::LinQuadSolverInstance, ::MOI.SimplexIterations) = true

# struct BarrierIterations <: AbstractSolverInstanceAttribute end
MOI.getattribute(m::LinQuadSolverInstance, ::MOI.BarrierIterations) = lqs_getbaritcnt(m.inner)
MOI.cangetattribute(m::LinQuadSolverInstance, ::MOI.BarrierIterations) = true

# struct NodeCount <: AbstractSolverInstanceAttribute end
MOI.getattribute(m::LinQuadSolverInstance, ::MOI.NodeCount) = lqs_getnodecnt(m.inner)
MOI.cangetattribute(m::LinQuadSolverInstance, ::MOI.NodeCount) = true

# struct RawSolver <: AbstractSolverInstanceAttribute end
MOI.getattribute(m::LinQuadSolverInstance, ::MOI.RawSolver) = m.inner
MOI.cangetattribute(m::LinQuadSolverInstance, ::MOI.RawSolver) = true