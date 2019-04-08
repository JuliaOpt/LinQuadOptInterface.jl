
mutable struct DataLinQuadOptimizer <: GenericLinQuadOptimizer
    LQOI.@LinQuadOptimizerBase
    params::Dict{String,Any}
    solve_function::Function
    DataLinQuadOptimizer(::Nothing) = new()
end

error_solve_function(instance) = error("No solve function set")

MOI.supports_constraint(m::LinQuadOptimizer, ::Type{F}, 
    ::Type{MOI.Interval}) where F <: MOI.AbstractFunction = false

function DataLinQuadOptimizer(
    solve_function::Function = error_solve_function; kwargs...)

    instance = DataLinQuadOptimizer(nothing)
    instance.params = Dict{String,Any}()
    instance.solve_function = solve_function

    MOI.empty!(instance)
    for (name, value) in kwargs
        instance.params[string(name)] = value
        setparam!(instance.inner, string(name), value)
    end
    return instance
end

mutable struct DataLinQuadOptimizerSolution
    termination_status::MOI.TerminationStatusCode

    primal_status::MOI.ResultStatusCode

    dual_status::MOI.ResultStatusCode

    variable_primal_solution::Vector{Float64}
    variable_dual_solution::Vector{Float64}

    constraint_primal_solution::Vector{Float64}
    constraint_dual_solution::Vector{Float64}

    quadratic_primal_solution::Vector{Float64}
    quadratic_dual_solution::Vector{Float64}

    objective_bound::Float64
    relative_mip_gap::Float64
    iteration_count::Int
    barrier_iterations::Int
    node_count::Int
end

get_parameters(instance::DataLinQuadOptimizer) = instance.params
get_optimization_sense(instance::DataLinQuadOptimizer) = instance.inner.sense
get_objective_coefficients(instance::DataLinQuadOptimizer) = instance.inner.c
get_objective_quadratic_coefficients(instance::DataLinQuadOptimizer) = instance.inner.Qobj
get_constraint_quadratic_coefficients(instance::DataLinQuadOptimizer) = instance.inner.Qcon
get_constraint_coefficients(instance::DataLinQuadOptimizer) = instance.inner.A
get_constraint_constant(instance::DataLinQuadOptimizer) = instance.inner.b
# get_constraint_constant_range_ub(instance::DataLinQuadOptimizer) = instance.inner.range
get_constraint_type(instance::DataLinQuadOptimizer) = instance.inner.contype
get_variable_lower_bound(instance::DataLinQuadOptimizer) = instance.inner.lb
get_variable_upper_bound(instance::DataLinQuadOptimizer) = instance.inner.ub
get_variable_type(instance::DataLinQuadOptimizer) = instance.inner.vartype
get_sos_constraints(instance::DataLinQuadOptimizer) = instance.inner.sos

function load_solution(instance::DataLinQuadOptimizer,
    solution::DataLinQuadOptimizerSolution)
    for field in fieldnames(DataLinQuadOptimizerSolution)
        setfield!(instance.inner, field, getfield(solution, field))
    end
    return
end

"""
    generic_solve(instance::DataLinQuadOptimizer)
"""
function generic_solve(instance::DataLinQuadOptimizer)
    solution = instance.solve_function(instance)::DataLinQuadOptimizerSolution
    load_solution(instance, solution)
    return
end
