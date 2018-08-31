function has_quadratic(model::LinQuadOptimizer)
    return model.obj_type == QuadraticObjective ||
        length(cmap(model).q_less_than) > 0 ||
        length(cmap(model).q_greater_than) > 0 ||
        length(cmap(model).q_equal_to) > 0
end

#=
    Optimize the model
=#

function MOI.optimize!(model::LinQuadOptimizer)
    # reset storage
    fill!(model.variable_primal_solution, NaN)
    fill!(model.variable_dual_solution, NaN)
    fill!(model.constraint_primal_solution, NaN)
    fill!(model.constraint_dual_solution, NaN)
    model.primal_status = MOI.UnknownResultStatus
    model.dual_status   = MOI.UnknownResultStatus
    model.primal_result_count = 0
    model.dual_result_count = 0

    start_time = time()
    if has_integer(model)
        solve_mip_problem!(model)
    elseif has_quadratic(model)
        solve_quadratic_problem!(model)
    else
        solve_linear_problem!(model)
    end
    model.solvetime = time() - start_time

    # termination_status
    model.termination_status = get_termination_status(model)
    model.primal_status = get_primal_status(model)
    model.dual_status = get_dual_status(model)

    if model.primal_status in [MOI.FeasiblePoint, MOI.InfeasiblePoint]
        get_variable_primal_solution!(model, model.variable_primal_solution)
        get_linear_primal_solution!(model, model.constraint_primal_solution)
        if has_quadratic(model)
            get_quadratic_primal_solution!(model, model.qconstraint_primal_solution)
        end
        model.primal_result_count = 1
    elseif model.primal_status == MOI.InfeasibilityCertificate
        get_unbounded_ray!(model, model.variable_primal_solution)
        model.primal_result_count = 1
    end
    if model.dual_status in [MOI.FeasiblePoint, MOI.InfeasiblePoint]
        get_variable_dual_solution!(model, model.variable_dual_solution)
        get_linear_dual_solution!(model, model.constraint_dual_solution)
        if has_quadratic(model)
            get_quadratic_dual_solution!(model, model.qconstraint_dual_solution)
        end
        model.dual_result_count = 1
    elseif model.dual_status == MOI.InfeasibilityCertificate
        get_farkas_dual!(model, model.constraint_dual_solution)
        get_farkas_dual_bounds!(model, model.variable_dual_solution)
        model.dual_result_count = 1
    end

    if MOI.get(model, MOI.ObjectiveSense()) == MOI.MaxSense
        model.constraint_dual_solution *= -1
        model.variable_dual_solution *= -1
    end
end


#=
    Result Count
=#
MOI.canget(::LinQuadOptimizer, ::MOI.ResultCount) = true
function MOI.get(model::LinQuadOptimizer, ::MOI.ResultCount)
    return max(model.primal_result_count, model.dual_result_count)
end

#=
    Termination status
=#
MOI.canget(::LinQuadOptimizer, ::MOI.TerminationStatus) = true
function MOI.get(model::LinQuadOptimizer, ::MOI.TerminationStatus)
    return model.termination_status
end

#=
    Primal status
=#
function MOI.canget(model::LinQuadOptimizer,  status::MOI.PrimalStatus)
    if status.N != 1
        error("Multiple primal statuses are not supported.")
    end
    return model.primal_result_count >= status.N
end
function MOI.get(model::LinQuadOptimizer, ::MOI.PrimalStatus)
    return model.primal_status
end

#=
    Dual status
=#
function MOI.canget(model::LinQuadOptimizer,  status::MOI.DualStatus)
    if status.N != 1
        error("Multiple dual statuses are not supported.")
    end
    return model.dual_result_count >= status.N
end
function MOI.get(model::LinQuadOptimizer, ::MOI.DualStatus)
    return model.dual_status
end

#=
    Objective Value
=#
function MOI.canget(::LinQuadOptimizer, attr::MOI.ObjectiveValue)
    return attr.resultindex == 1
end
function MOI.get(model::LinQuadOptimizer, attr::MOI.ObjectiveValue)
    if attr.resultindex == 1
        # Note: we add m.objective_constant here to account for any constant
        # term which was not passed to the solver itself (and which therefore
        # would not be accounted for in `get_objective_value(m)`. We do *not*
        # call `get_constant_objective(m)` because that would also pull any
        # constants which were passed to the solver, resulting those constants
        # being counted twice. This confusion will be alleviated when all LQOI
        # solvers implement `get_constant_objective()` and
        # `set_constant_objective!()` by actually passing the relevant constants
        # to the solvers, at which point we can just get rid of
        # m.objective_constant entirely.
        return get_objective_value(model) + model.objective_constant
    else
        error("Unable to access multiple objective values")
    end
end

#=
    Variable Primal solution
=#
MOI.canget(::LinQuadOptimizer, ::MOI.VariablePrimal, ::Type{VarInd}) = true
function MOI.get(model::LinQuadOptimizer, ::MOI.VariablePrimal, index::VarInd)
    column = get_column(model, index)
    return model.variable_primal_solution[column]
end

MOI.canget(::LinQuadOptimizer, ::MOI.VariablePrimal, ::Type{Vector{VarInd}}) = true
function MOI.get(model::LinQuadOptimizer, ::MOI.VariablePrimal, indices::Vector{VarInd})
    MOI.get.(Ref(model), Ref(MOI.VariablePrimal()), indices)
end

#=
    Variable Dual solution
=#
"""
    is_binding(set, value::Float64)

Return true if `value` is an extreme point of the set `set`.
"""
is_binding(set::LE, value::Float64) = isapprox(set.upper, value)
is_binding(set::GE, value::Float64) = isapprox(set.lower, value)
is_binding(set::EQ, value::Float64) = isapprox(set.value, value)
is_binding(set::IV, value::Float64) = isapprox(set.lower, value) || isapprox(set.upper, value)

function MOI.canget(::LinQuadOptimizer, ::MOI.ConstraintDual, ::Type{SVCI{S}}) where S <: LinSets
    return true
end

function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintDual, index::SVCI{<: LinSets})
    column = get_column(model, model[index])
    # the variable reduced cost is only the constraint dual if the bound is active.
    # or it might be a dual ray
    if model.dual_status == MOI.InfeasibilityCertificate
        return model.variable_dual_solution[column]
    else
        set = MOI.get(model, MOI.ConstraintSet(), index)
        primal_value = model.variable_primal_solution[column]
        if is_binding(set, primal_value)
            return model.variable_dual_solution[column]
        else
            return 0.0
        end
    end
 end

function MOI.canget(s::LinQuadOptimizer, ::MOI.ConstraintDual,
                    ::Type{VVCI{S}}) where S <: Union{MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives}
    return true
end
function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintDual, index::VVCI{<: Union{MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives}})
    return [model.constraint_dual_solution[row] for row in model[index]]
end

#=
    Variable Bound Primal solution
=#

MOI.canget(::LinQuadOptimizer, ::MOI.ConstraintPrimal, ::Type{SVCI{S}}) where S <: LinSets = true

function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintPrimal, index::SVCI{<: LinSets})
    column = get_column(model, model[index])
    return model.variable_primal_solution[column]
end

function MOI.canget(::LinQuadOptimizer, ::MOI.ConstraintPrimal,
                   ::Type{VVCI{S}}) where S <: Union{MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives}
    return true
end

function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintPrimal,
                 index::VVCI{<: Union{MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives}})
    return [model.constraint_primal_solution[row] for row in model[index]]
end

#=
    Constraint Primal solution
=#

MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintPrimal, ::Type{LCI{S}}) where S <: LinSets = true

function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintPrimal, index::LCI{<: LinSets})
    row = model[index]
    return model.constraint_primal_solution[row] + model.constraint_constant[row]
end

MOI.canget(::LinQuadOptimizer, ::MOI.ConstraintPrimal, ::Type{QCI{S}}) where S <: LinSets = true

function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintPrimal, index::QCI{<: LinSets})
    row = model[index]
    return model.qconstraint_primal_solution[row]
end

MOI.canget(::LinQuadOptimizer, ::MOI.ConstraintPrimal, ::Type{VLCI{S}}) where S <: Union{MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives} = true

# vector valued constraint duals
function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintPrimal, index::VLCI{<: Union{MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives}})
    row = model[index]
    return model.constraint_primal_solution[row] + model.constraint_constant[row]
end

#=
    Constraint Dual solution
=#

__assert_dual_sense__(::LCI{LE}, dual) = @assert dual <= 0.0
__assert_dual_sense__(::LCI{GE}, dual) = @assert dual >= 0.0
__assert_dual_sense__(::LCI{IV}, dual) = nothing
__assert_dual_sense__(::LCI{EQ}, dual) = nothing

MOI.canget(::LinQuadOptimizer, ::MOI.ConstraintDual, ::Type{LCI{S}}) where S <: LinSets = true

function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintDual, index::LCI{<: LinSets})
    row = model[index]
    dual = model.constraint_dual_solution[row]
    __assert_dual_sense__(index, dual)
    return dual
end

MOI.canget(model::LinQuadOptimizer, ::MOI.ConstraintDual, ::Type{QCI{S}}) where S <: LinSets = true

function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintDual, index::QCI{<: LinSets})
    row = model[index]
    return model.qconstraint_dual_solution[row]
end


# vector valued constraint duals
MOI.canget(::LinQuadOptimizer, ::MOI.ConstraintDual, ::Type{VLCI{S}}) where S <: Union{MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives} = true

function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintDual, index::VLCI{<: Union{MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives}})
    rows = model[index]
    return model.constraint_dual_solution[rows]
end

#=
    Solution Attributes
=#

MOI.supports(::LinQuadOptimizer, ::MOI.ObjectiveBound) = true
MOI.canget(::LinQuadOptimizer, ::MOI.ObjectiveBound) = true
function MOI.get(model::LinQuadOptimizer, ::MOI.ObjectiveBound)
    return get_objective_bound(model)
end

MOI.supports(::LinQuadOptimizer, ::MOI.RelativeGap) = true
MOI.canget(::LinQuadOptimizer, ::MOI.RelativeGap) = true
function MOI.get(model::LinQuadOptimizer, ::MOI.RelativeGap)
    return get_relative_mip_gap(model)
end

MOI.supports(::LinQuadOptimizer, ::MOI.SolveTime) = true
MOI.canget(::LinQuadOptimizer, ::MOI.SolveTime) = true
function MOI.get(model::LinQuadOptimizer, ::MOI.SolveTime)
    return model.solvetime
end

MOI.supports(::LinQuadOptimizer, ::MOI.SimplexIterations) = true
MOI.canget(::LinQuadOptimizer, ::MOI.SimplexIterations) = true
function MOI.get(model::LinQuadOptimizer, ::MOI.SimplexIterations)
    return get_iteration_count(model)
end

MOI.supports(::LinQuadOptimizer, ::MOI.BarrierIterations) = true
MOI.canget(::LinQuadOptimizer, ::MOI.BarrierIterations) = true
function MOI.get(model::LinQuadOptimizer, ::MOI.BarrierIterations)
    return get_barrier_iterations(model)
end

MOI.supports(::LinQuadOptimizer, ::MOI.NodeCount) = true
MOI.canget(::LinQuadOptimizer, ::MOI.NodeCount) = true
function MOI.get(model::LinQuadOptimizer, ::MOI.NodeCount)
    return get_node_count(model)
end

MOI.supports(::LinQuadOptimizer, ::MOI.RawSolver) = true
MOI.canget(::LinQuadOptimizer, ::MOI.RawSolver) = true
function MOI.get(model::LinQuadOptimizer, ::MOI.RawSolver)
    return model
end
