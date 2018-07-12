function hasquadratic(m::LinQuadOptimizer)
    (m.obj_type == QuadraticObjective) || (length(cmap(m).q_less_than) + length(cmap(m).q_greater_than) + length(cmap(m).q_equal_to) > 0)
end

#=
    Optimize the model
=#

function MOI.optimize!(m::LinQuadOptimizer)
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
        solve_mip_problem!(m)
    elseif hasquadratic(m)
        solve_quadratic_problem!(m)
    else
        solve_linear_problem!(m)
    end
    m.solvetime = time() - t

    # termination_status
    m.termination_status = get_termination_status(m)
    m.primal_status = get_primal_status(m)
    m.dual_status = get_dual_status(m)

    if m.primal_status in [MOI.FeasiblePoint, MOI.InfeasiblePoint]
        # primal solution exists
        get_variable_primal_solution!(m, m.variable_primal_solution)
        get_linear_primal_solution!(m, m.constraint_primal_solution)
        if hasquadratic(m)
            get_quadratic_primal_solution!(m, m.qconstraint_primal_solution)
        end
        m.primal_result_count = 1
        # CPLEX can return infeasible points
    elseif m.primal_status == MOI.InfeasibilityCertificate
        get_unbounded_ray!(m, m.variable_primal_solution)
        m.primal_result_count = 1
    end
    if m.dual_status in [MOI.FeasiblePoint, MOI.InfeasiblePoint]
        # dual solution exists
        get_variable_dual_solution!(m, m.variable_dual_solution)
        get_linear_dual_solution!(m, m.constraint_dual_solution)
        if hasquadratic(m)
            get_quadratic_dual_solution!(m, m.qconstraint_dual_solution)
        end
        m.dual_result_count = 1
        # dual solution may not be feasible
    elseif m.dual_status == MOI.InfeasibilityCertificate
        get_farkas_dual!(m, m.constraint_dual_solution)
        m.dual_result_count = 1
    end

    #=
        CPLEX has the dual convention that the sign of the dual depends on the
        optimization sense. This isn't the same as the MOI convention so we need
        to correct that.
    =#
    # TODO
    if MOI.get(m, MOI.ObjectiveSense()) == MOI.MaxSense
        m.constraint_dual_solution *= -1
        m.variable_dual_solution *= -1
    end
end


#=
    Result Count
=#
function MOI.get(m::LinQuadOptimizer, ::MOI.ResultCount)
    max(m.primal_result_count, m.dual_result_count)
end
MOI.canget(m::LinQuadOptimizer, ::MOI.ResultCount) = true

#=
    Termination status
=#

function MOI.get(m::LinQuadOptimizer, ::MOI.TerminationStatus)
    m.termination_status
end
MOI.canget(m::LinQuadOptimizer, ::MOI.TerminationStatus) = true

#=
    Primal status
=#

function MOI.get(m::LinQuadOptimizer, p::MOI.PrimalStatus)
    m.primal_status
end
function MOI.canget(m::LinQuadOptimizer, p::MOI.PrimalStatus)
    m.primal_result_count >= p.N
end

#=
    Dual status
=#

function MOI.get(m::LinQuadOptimizer, d::MOI.DualStatus)
    m.dual_status
end
function MOI.canget(m::LinQuadOptimizer, d::MOI.DualStatus)
    m.dual_result_count >= d.N
end

#=
    Objective Value
=#


function MOI.get(m::LinQuadOptimizer, attr::MOI.ObjectiveValue)
    if attr.resultindex == 1
        get_objective_value(m) + m.objective_constant
    else
        error("Unable to access multiple objective values")
    end
end
function MOI.canget(m::LinQuadOptimizer, attr::MOI.ObjectiveValue)
    if attr.resultindex == 1
        return true
    else
        return false
    end
end

#=
    Variable Primal solution
=#


function MOI.get(m::LinQuadOptimizer, ::MOI.VariablePrimal, v::VarInd)
    col = m.variable_mapping[v]
    return m.variable_primal_solution[col]
end
MOI.canget(m::LinQuadOptimizer, ::MOI.VariablePrimal, ::Type{VarInd}) = true

function MOI.get(m::LinQuadOptimizer, ::MOI.VariablePrimal, v::Vector{VarInd})
    MOI.get.(Ref(m), MOI.VariablePrimal(), v)
end
MOI.canget(m::LinQuadOptimizer, ::MOI.VariablePrimal, ::Type{Vector{VarInd}}) = true

#=
    Variable Dual solution
=#
isbinding(set::LE, value::Float64) = isapprox(set.upper, value)
isbinding(set::GE, value::Float64) = isapprox(set.lower, value)
isbinding(set::EQ, value::Float64) = isapprox(set.value, value)
isbinding(set::IV, value::Float64) = isapprox(set.lower, value) || isapprox(set.upper, value)
function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintDual, c::SVCI{<: LinSets})
    vref = m[c]
    col = m.variable_mapping[vref]
    # the variable reduced cost is only the constriant dual if the bound is active.
    set = MOI.get(m, MOI.ConstraintSet(), c)
    solval = m.variable_primal_solution[col]
    if isbinding(set, solval)
        return m.variable_dual_solution[col]
    else
        return 0.0
    end
end
MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintDual, ::Type{SVCI{S}}) where S <: LinSets = true

function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintDual, c::VVCI{<: Union{MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives}})
    rowrefs = m[c]
    out = Float64[]
    sizehint!(out, length(rowrefs))
    for ref in rowrefs
        push!(out, m.constraint_dual_solution[ref])
    end
    return out
end
MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintDual, ::Type{VVCI{S}}) where S <: Union{MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives} = true

#=
    Variable Bound Primal solution
=#

function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintPrimal, c::SVCI{<: LinSets})
    vref = m[c]
    col = m.variable_mapping[vref]
    return m.variable_primal_solution[col]
end
MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintPrimal, ::Type{SVCI{S}}) where S <: LinSets = true

function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintPrimal, c::VVCI{<: Union{MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives}})
    rowrefs = m[c]
    out = Float64[]
    sizehint!(out, length(rowrefs))
    for ref in rowrefs
        push!(out, m.constraint_primal_solution[ref])
    end
    return out
end
MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintPrimal, ::Type{VVCI{S}}) where S <: Union{MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives} = true

#=
    Constraint Primal solution
=#

function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintPrimal, c::LCI{<: LinSets})
    row = m[c]
    return m.constraint_primal_solution[row]+m.constraint_constant[row]
end
MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintPrimal, ::Type{LCI{S}}) where S <: LinSets = true

function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintPrimal, c::QCI{<: LinSets})
    row = m[c]
    return m.qconstraint_primal_solution[row]#+m.qconstraint_constant[row]
end
MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintPrimal, ::Type{QCI{S}}) where S <: LinSets = true

# vector valued constraint duals
function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintPrimal, c::VLCI{<: Union{MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives}})
    row = m[c]
    return m.constraint_primal_solution[row]+m.constraint_constant[row]
end
MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintPrimal, ::Type{VLCI{S}}) where S <: Union{MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives} = true

#=
    Constraint Dual solution
=#

_checkdualsense(::LCI{LE}, dual) = dual <= 0.0
_checkdualsense(::LCI{GE}, dual) = dual >= 0.0
_checkdualsense(::LCI{IV}, dual) = true
_checkdualsense(::LCI{EQ}, dual) = true

function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintDual, c::LCI{<: LinSets})
    row = m[c]
    dual = m.constraint_dual_solution[row]
    @assert _checkdualsense(c, dual)
    return dual
end
MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintDual, ::Type{LCI{S}}) where S <: LinSets = true

function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintDual, c::QCI{<: LinSets})
    row = m[c]
    return m.qconstraint_dual_solution[row]#+m.qconstraint_constant[row]
end
MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintDual, ::Type{QCI{S}}) where S <: LinSets = true


# vector valued constraint duals
MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintDual, c::VLCI{<: Union{MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives}}) = m.constraint_dual_solution[m[c]]
MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintDual, ::Type{VLCI{S}}) where S <: Union{MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives} = true

#=
    Solution Attributes
=#

# struct ObjectiveBound <: MOI.AbstractOptimizerAttribute end
MOI.get(m::LinQuadOptimizer, ::MOI.ObjectiveBound) = get_objective_bound(m)
MOI.canget(m::LinQuadOptimizer, ::MOI.ObjectiveBound) = true

# struct RelativeGap <: MOI.AbstractOptimizerAttribute  end
MOI.get(m::LinQuadOptimizer, ::MOI.RelativeGap) = get_relative_mip_gap(m)
MOI.canget(m::LinQuadOptimizer, ::MOI.RelativeGap) = true

# struct SolveTime <: MOI.AbstractOptimizerAttribute end
MOI.get(m::LinQuadOptimizer, ::MOI.SolveTime) = m.solvetime
MOI.canget(m::LinQuadOptimizer, ::MOI.SolveTime) = true

# struct SimplexIterations <: MOI.AbstractOptimizerAttribute end
MOI.get(m::LinQuadOptimizer, ::MOI.SimplexIterations) = get_iteration_count(m)
MOI.canget(m::LinQuadOptimizer, ::MOI.SimplexIterations) = true

# struct BarrierIterations <: MOI.AbstractOptimizerAttribute end
MOI.get(m::LinQuadOptimizer, ::MOI.BarrierIterations) = get_barrier_iterations(m)
MOI.canget(m::LinQuadOptimizer, ::MOI.BarrierIterations) = true

# struct NodeCount <: MOI.AbstractOptimizerAttribute end
MOI.get(m::LinQuadOptimizer, ::MOI.NodeCount) = get_node_count(m)
MOI.canget(m::LinQuadOptimizer, ::MOI.NodeCount) = true

# struct RawSolver <: MOI.AbstractOptimizerAttribute end
MOI.get(m::LinQuadOptimizer, ::MOI.RawSolver) = m
MOI.canget(m::LinQuadOptimizer, ::MOI.RawSolver) = true
