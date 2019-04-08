
mutable struct MockLinQuadOptimizer <: GenericLinQuadOptimizer
    LQOI.@LinQuadOptimizerBase
    params::Dict{String,Any}

    termination_status_stored::Vector{MOI.TerminationStatusCode}
    primal_status_stored::Vector{MOI.ResultStatusCode}
    dual_status_stored::Vector{MOI.ResultStatusCode}

    variable_primal_solution_stored::Vector{Vector{Float64}}
    variable_dual_solution_stored::Vector{Vector{Float64}}

    constraint_primal_solution_stored::Vector{Vector{Float64}}
    constraint_dual_solution_stored::Vector{Vector{Float64}}

    quadratic_primal_solution_stored::Vector{Vector{Float64}}
    quadratic_dual_solution_stored::Vector{Vector{Float64}}

    MockLinQuadOptimizer(::Nothing) = new()
end

function MockLinQuadOptimizer(;kwargs...)

    instance = MockLinQuadOptimizer(nothing)
    instance.params = Dict{String,Any}()

    instance.termination_status_stored = MOI.TerminationStatusCode[]
    instance.primal_status_stored = MOI.ResultStatusCode[]
    instance.dual_status_stored = MOI.ResultStatusCode[]

    instance.variable_primal_solution_stored = Vector{Float64}[]
    instance.variable_dual_solution_stored = Vector{Float64}[]

    instance.constraint_primal_solution_stored = Vector{Float64}[]
    instance.constraint_dual_solution_stored = Vector{Float64}[]

    instance.quadratic_primal_solution_stored = Vector{Float64}[]
    instance.quadratic_dual_solution_stored = Vector{Float64}[]

    MOI.empty!(instance)
    for (name, value) in kwargs
        instance.params[string(name)] = value
        setparam!(instance.inner, string(name), value)
    end
    return instance
end

"""
    unload(from::Vector, to, warn = true)

Helper function to remove the first element of a vector and return it.
If the vector is empty data in `default` is returned instead.
Used in `fakesolve`.
"""
function unload(from, default, warn = true)
    if !isempty(from)
        out = from[1]
        Compat.popfirst!(from)
        return out
    else
        if warn
            Compat.@warn("No data in the input vector, returning default.")
        end
        return default
    end
end

"""
    fakesolve(instance::MockLinQuadOptimizer)

Set solutions upon solve calls.
Data held in `MockLinQuadOptimizer` stored data is passed to the low-level emulator `GenericLinQuadModel`.
"""
function generic_solve(instance::MockLinQuadOptimizer)

    instance.inner.termination_status = unload(instance.termination_status_stored, instance.inner.termination_status)
    instance.inner.primal_status = unload(instance.primal_status_stored, instance.inner.primal_status)
    instance.inner.dual_status = unload(instance.dual_status_stored, instance.inner.dual_status)

    instance.inner.variable_primal_solution = unload(instance.variable_primal_solution_stored, instance.inner.variable_primal_solution)
    instance.inner.variable_dual_solution = unload(instance.variable_dual_solution_stored, instance.inner.variable_dual_solution)

    instance.inner.constraint_primal_solution = unload(instance.constraint_primal_solution_stored, instance.inner.constraint_primal_solution)
    instance.inner.constraint_dual_solution = unload(instance.constraint_dual_solution_stored, instance.inner.constraint_dual_solution)

    instance.inner.quadratic_primal_solution = unload(instance.quadratic_primal_solution_stored, instance.inner.quadratic_primal_solution)
    instance.inner.quadratic_dual_solution = unload(instance.quadratic_dual_solution_stored, instance.inner.quadratic_dual_solution)

    try
        instance.inner.objective_bound = get_objval(instance.inner)
    catch
        instance.inner.objective_bound = NaN
    end
    return
end

function set_variable_primal_solution!(instance::MockLinQuadOptimizer, input)
    push!(instance.variable_primal_solution_stored, input)
end
function set_variable_dual_solution!(instance::MockLinQuadOptimizer, input)
    push!(instance.variable_dual_solution_stored, input)
end
function set_constraint_primal_solution!(instance::MockLinQuadOptimizer, input)
    push!(instance.constraint_primal_solution_stored, input)
end
function set_constraint_dual_solution!(instance::MockLinQuadOptimizer, input)
    push!(instance.constraint_dual_solution_stored, input)
end
function set_quadratic_dual_solution!(instance::MockLinQuadOptimizer, input)
    push!(instance.quadratic_dual_solution_stored, input)
end
function set_quadratic_primal_solution!(instance::MockLinQuadOptimizer, input)
    push!(instance.quadratic_primal_solution_stored, input)
end
function set_termination_status!(instance::MockLinQuadOptimizer, input)
    push!(instance.termination_status_stored, input)
end
function set_primal_status!(instance::MockLinQuadOptimizer, input)
    push!(instance.primal_status_stored, input)
end
function set_dual_status!(instance::MockLinQuadOptimizer, input)
    push!(instance.dual_status_stored, input)
end

function set_solution!(instance::MockLinQuadOptimizer;
    variable_primal = Float64[],
    variable_dual = Float64[],
    constraint_primal = Float64[],
    constraint_dual = Float64[],
    quadratic_primal = Float64[],
    quadratic_dual = Float64[],
    termination_status = MOI.OPTIMAL,
    primal_status = MOI.FEASIBLE_POINT,
    dual_status = MOI.FEASIBLE_POINT
    )

    set_variable_primal_solution!(instance, variable_primal)
    set_variable_dual_solution!(instance, variable_dual)
    set_constraint_primal_solution!(instance, constraint_primal)
    set_constraint_dual_solution!(instance, constraint_dual)
    set_quadratic_primal_solution!(instance, quadratic_primal)
    set_quadratic_dual_solution!(instance, quadratic_dual)
    set_termination_status!(instance, termination_status)
    set_primal_status!(instance, primal_status)
    set_dual_status!(instance, dual_status)

    return
end
