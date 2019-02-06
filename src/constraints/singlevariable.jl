# Notes to begin.
#
# For all SingleVariable-in-X constraints, we use the same index as the variable
# index. This implies that we cannot add two constraints of the same type. In
# addition, some combinations of constraints are forbidden. Unfortunately, this
# means that SingleVariable-in-Sets are not returned ordered  by creation with
# ListOfConstraintIndices, they are ordered by variable creation.

function throw_conflicting_variable_type(set, variable_type)
    error("Cannot add $(set) constraint because variable is $(variable_type).")
end

function throw_conflicting_upper_bound(set)
    error("Cannot add $(set) constraint because an upper bound already exists.")
end

function throw_conflicting_lower_bound(set)
    error("Cannot add $(set) constraint because a lower bound already exists.")
end

"""
    has_constraint_filter(cache::VariableCache, set::S) where {S}

Return `true` if the variable associated with the variable cache `cache` has a
constraint of type `MOI.ConstraintIndex{MOI.SingleVariable, S}`.
"""
has_constraint_filter(::VariableCache, arg) = false

"""
    variable_index(c_index::SVCI)

Convert a `ConstraintIndex{SingleVariable, S}` to a `VariableIndex`.
"""
function variable_index(c_index::SVCI)
    return MOI.VariableIndex(c_index.value)
end

variable_index(index::MOI.SingleVariable) = index.variable

"""
    single_variable(c_index::SVCI)

Convert a `ConstraintIndex{SingleVariable, S}` to a `SingleVariable`.
"""
function single_variable(c_index::SVCI)
    return MOI.SingleVariable(variable_index(c_index))
end

"""
    constraint_index(variable::VariableIndex, ::S)

Convert a `VariableIndex` to a `ConstraintIndex{SingleVariable, S}`.
"""
function constraint_index(variable::MOI.VariableIndex, ::S) where {S}
    return MOI.ConstraintIndex{MOI.SingleVariable, S}(variable.value)
end

"""
    constraint_index(variable::SingleVariable, ::S)

Convert a `SingleVariable` to a `ConstraintIndex{SingleVariable, S}`.
"""
function constraint_index(variable::MOI.SingleVariable, set::S) where {S}
    return constraint_index(variable.variable, set)
end

function variable_cache(model::LinQuadOptimizer, index::MOI.VariableIndex)
    return model.variable_cache[index]
end

function variable_cache(model::LinQuadOptimizer, index::SVCI)
    return variable_cache(model, variable_index(index))
end

function variable_cache(model::LinQuadOptimizer, index::MOI.SingleVariable)
    return variable_cache(model, variable_index(index))
end

function MOI.is_valid(model::LinQuadOptimizer, index::SVCI{S}) where {S}
    cache = variable_cache(model, index)
    return has_constraint_filter(cache, S)
end

function MOI.get(::LinQuadOptimizer, ::MOI.ConstraintFunction, index::SVCI)
    return single_variable(index)
end

# TODO(odow): we could cache this. It seems very inefficient.
function MOI.get(model::LinQuadOptimizer, ::MOI.NumberOfConstraints{MOI.SingleVariable, S}) where {S}
    num = 0
    for (index, cache) in model.variable_cache
        if has_constraint_filter(cache, S)
            num += 1
        end
    end
    return num
end

# Get the list of constraint indices. We order them by variable creation time,
# i.e., the order of variables in model.variable_references, which should be
# their column ordering.
function MOI.get(model::LinQuadOptimizer, ::MOI.ListOfConstraintIndices{MOI.SingleVariable, S}) where {S}
    indices = MOI.ConstraintIndex{MOI.SingleVariable, S}[]
    for variable in model.variable_references
        cache = variable_cache(model, variable)
        if has_constraint_filter(cache, S)
            push!(indices, MOI.ConstraintIndex{MOI.SingleVariable, S}(variable.value))
        end
    end
    return indices
end

###
### SingleVariable in LessThan, GreaterThan, Interval, EqualTo
###

function set_cache(cache::VariableCache, set::MOI.LessThan)
    cache.upper = set.upper
    cache.bound_type = cache.bound_type == ONLY_GREATER_THAN ? LESS_AND_GREATER_THAN : ONLY_LESS_THAN
    return
end

function unset_cache(cache::VariableCache, ::Type{<:MOI.LessThan})
    cache.upper = Inf
    cache.bound_type = cache.bound_type == LESS_AND_GREATER_THAN ? ONLY_GREATER_THAN : FREE
    return
end

function set_cache(cache::VariableCache, set::MOI.GreaterThan)
    cache.lower = set.lower
    cache.bound_type = cache.bound_type == ONLY_LESS_THAN ? LESS_AND_GREATER_THAN : ONLY_GREATER_THAN
    return
end

function unset_cache(cache::VariableCache, ::Type{<:MOI.GreaterThan})
    cache.lower = -Inf
    cache.bound_type = cache.bound_type == LESS_AND_GREATER_THAN ? ONLY_LESS_THAN : FREE
    return
end

function set_cache(cache::VariableCache, set::MOI.Interval)
    cache.lower = set.lower
    cache.upper = set.upper
    cache.bound_type = INTERVAL
    return
end

function set_cache(cache::VariableCache, set::MOI.EqualTo)
    cache.lower = set.value
    cache.upper = set.value
    cache.bound_type = EQUAL_TO
    return
end

function unset_cache(cache::VariableCache, ::Type{<:Union{MOI.Interval, MOI.EqualTo}})
    cache.lower = -Inf
    cache.upper = Inf
    cache.bound_type = FREE
    return
end

function MOI.add_constraint(model::LinQuadOptimizer, variable::SinVar, set::S) where S <: LinSets
    __assert_supported_constraint__(model, SinVar, S)
    cache = variable_cache(model, variable)
    check_conflicting_constraint(cache, set)
    set_variable_bound(model, variable, set)
    set_cache(cache, set)
    return constraint_index(variable, set)
end

function MOI.add_constraints(model::LinQuadOptimizer, variables::Vector{SinVar},
                             sets::Vector{S}) where S <: LinSets
    __assert_supported_constraint__(model, SinVar, S)
    for (variable, set) in zip(variables, sets)
        check_conflicting_constraint(variable_cache(model, variable), set)
    end
    set_variable_bounds(model, variables, sets)
    indices = SVCI{S}[]
    for (variable, set) in zip(variables, sets)
        cache = variable_cache(model, variable)
        set_cache(cache, set)
        push!(indices, constraint_index(variable, set))
    end
    return indices
end

function MOI.delete(model::LinQuadOptimizer, index::SVCI{S}) where {S <: LinSets}
    __assert_valid__(model, index)
    cache = variable_cache(model, index)
    unset_cache(cache, S)
    set_variable_bound(model, single_variable(index), IV(cache.lower, cache.upper))
    return
end

function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintSet, index::SVCI{LE})
    cache = variable_cache(model, index)
    return MOI.LessThan{Float64}(cache.upper)
end

function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintSet, index::SVCI{GE})
    cache = variable_cache(model, index)
    return MOI.GreaterThan{Float64}(cache.lower)
end

function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintSet, index::SVCI{EQ})
    cache = variable_cache(model, index)
    return MOI.EqualTo{Float64}(cache.lower)
end

function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintSet, index::SVCI{IV})
    cache = variable_cache(model, index)
    return MOI.Interval{Float64}(cache.lower, cache.upper)
end

MOI.supports(::LinQuadOptimizer, ::MOI.ConstraintSet, ::Type{SVCI{<:LinSets}}) = true
function MOI.set(model::LinQuadOptimizer, ::MOI.ConstraintSet, index::SVCI{S}, new_set::S) where S <: LinSets
    cache = variable_cache(model, index)
    set_variable_bound(model, single_variable(index), new_set)
    set_cache(cache, new_set)
    return
end

function has_constraint_filter(cache::VariableCache, ::Type{<:MOI.LessThan})
    return cache.bound_type == ONLY_LESS_THAN || cache.bound_type == LESS_AND_GREATER_THAN
end

function has_constraint_filter(cache::VariableCache, ::Type{<:MOI.GreaterThan})
    return cache.bound_type == ONLY_GREATER_THAN || cache.bound_type == LESS_AND_GREATER_THAN
end

function has_constraint_filter(cache::VariableCache, ::Type{<:MOI.EqualTo})
    return cache.bound_type == EQUAL_TO
end

function has_constraint_filter(cache::VariableCache, ::Type{<:MOI.Interval})
    return cache.bound_type == INTERVAL
end

function check_conflicting_constraint(cache::VariableCache, ::MOI.LessThan)
    if cache.variable_type == SEMI_INTEGER || cache.variable_type == SEMI_CONTINUOUS
        throw_conflicting_variable_type(MOI.LessThan, cache.variable_type)
    elseif cache.upper != Inf
        throw_conflicting_upper_bound(MOI.LessThan)
    end
    return
end

function check_conflicting_constraint(cache::VariableCache, ::MOI.GreaterThan)
    if cache.variable_type == SEMI_INTEGER || cache.variable_type == SEMI_CONTINUOUS
        throw_conflicting_variable_type(MOI.GreaterThan, cache.variable_type)
    elseif cache.lower != -Inf
        throw_conflicting_lower_bound(MOI.GreaterThan)
    end
    return
end

function check_conflicting_constraint(cache::VariableCache, ::MOI.Interval)
    if cache.variable_type == SEMI_INTEGER || cache.variable_type == SEMI_CONTINUOUS || cache.variable_type == BINARY
        throw_conflicting_variable_type(MOI.Interval, cache.variable_type)
    elseif cache.lower != -Inf
        throw_conflicting_lower_bound(MOI.Interval)
    elseif cache.upper != Inf
        throw_conflicting_upper_bound(MOI.Interval)
    end
    return
end

function check_conflicting_constraint(cache::VariableCache, ::MOI.EqualTo)
    if cache.variable_type == SEMI_INTEGER || cache.variable_type == SEMI_CONTINUOUS
        throw_conflicting_variable_type(MOI.EqualTo, cache.variable_type)
    elseif cache.lower != -Inf
        throw_conflicting_lower_bound(MOI.EqualTo)
    elseif cache.upper != Inf
        throw_conflicting_upper_bound(MOI.EqualTo)
    end
    return
end

"""
    change_both_variable_bounds!(model::LinQuadOptimizer, columns::Vector{Int},
        lower_bounds::Vector{Float64}, upper_bounds::Vector{Float64})

Set the lower bound of column `columns[i]` to `lower_bounds[i]` and the upper
bound to `upper_bounds[i]`. Alternatively, the lower or upper bound can be left
blank by passing an array of length 0 instead.

Examples:
    change_both_variable_bounds!(model, [1, 2], [-0.5, 0.0], [1.0, 2.0])
    change_both_variable_bounds!(model, [1, 2], [-0.5, 0.0], Float64[])
    change_both_variable_bounds!(model, [1, 2], Float64[], [1.0, 2.0])
"""
function change_both_variable_bounds!(
        model::LinQuadOptimizer,
        columns::Vector{Int},
        lower_bounds::Vector{Float64},
        upper_bounds::Vector{Float64})
    if length(lower_bounds) > 0 && length(upper_bounds) > 0
        columns = vcat(columns, columns)
    end
    change_variable_bounds!(model,
        columns,
        vcat(lower_bounds, upper_bounds),
        vcat(
            fill(backend_type(model, Val{:Lowerbound}()), length(lower_bounds)),
            fill(backend_type(model, Val{:Upperbound}()), length(upper_bounds))
        )
    )
end

function set_variable_bound(model::LinQuadOptimizer, v::SinVar, set::LE)
    change_both_variable_bounds!(
        model, [get_column(model, v)], Float64[], [set.upper])
end

function set_variable_bound(model::LinQuadOptimizer, v::SinVar, set::GE)
    change_both_variable_bounds!(
        model, [get_column(model, v)], [set.lower], Float64[])
end

function set_variable_bound(model::LinQuadOptimizer, v::SinVar, set::EQ)
    change_both_variable_bounds!(
        model, [get_column(model, v)], [set.value], [set.value])
end

function set_variable_bound(model::LinQuadOptimizer, v::SinVar, set::IV)
    change_both_variable_bounds!(
        model, [get_column(model, v)], [set.lower], [set.upper])
end

function set_variable_bounds(
        model::LinQuadOptimizer, vars::Vector{SinVar}, sets::Vector{LE})
    change_both_variable_bounds!(model, get_column.(model, vars), Float64[],
        [set.upper for set in sets])
end

function set_variable_bounds(
        model::LinQuadOptimizer, vars::Vector{SinVar}, sets::Vector{GE})
    change_both_variable_bounds!(model, get_column.(model, vars),
        [set.lower for set in sets], Float64[])
end

function set_variable_bounds(
        model::LinQuadOptimizer, vars::Vector{SinVar}, sets::Vector{EQ})
    values = [set.value for set in sets]
    change_both_variable_bounds!(model, get_column.(model, vars), values, values)
end

function set_variable_bounds(
        model::LinQuadOptimizer, vars::Vector{SinVar}, sets::Vector{IV})
    change_both_variable_bounds!(model, get_column.(model, vars),
        [set.lower for set in sets], [set.upper for set in sets])
end

###
### SingleVariable in ZeroOne
###
### Conflicts with Integer, Semicontinuous, and Semiinteger constraints.
###
### Note: we cache existing bounds.

function MOI.add_constraint(model::LinQuadOptimizer, variable::SinVar, set::MOI.ZeroOne)
    __assert_supported_constraint__(model, SinVar, MOI.ZeroOne)
    cache = variable_cache(model, variable)
    check_conflicting_constraint(cache, set)
    cache.variable_type = BINARY
    cache.lower_old = cache.lower
    cache.upper_old = cache.upper
    cache.lower = max(0.0, cache.lower)
    cache.upper = min(1.0, cache.upper)
    set_variable_bound(model, variable, MOI.Interval(cache.lower, cache.upper))
    change_variable_types!(model, [cache.column], [backend_type(model, set)])
    make_problem_type_integer(model)
    return constraint_index(variable, set)
end

function MOI.delete(model::LinQuadOptimizer, index::SVCI{MOI.ZeroOne})
    __assert_valid__(model, index)
    cache = variable_cache(model, index)
    cache.variable_type = CONTINUOUS
    cache.lower = cache.lower_old
    cache.upper = cache.upper_old
    set_variable_bound(model, single_variable(index), MOI.Interval(cache.lower, cache.upper))
    change_variable_types!(
        model, [cache.column], [backend_type(model, Val{:Continuous}())])
    if !has_integer(model)
        make_problem_type_continuous(model)
    end
end

function MOI.get(::LinQuadOptimizer, ::MOI.ConstraintSet, ::SVCI{MOI.ZeroOne})
    return MOI.ZeroOne()
end

function has_constraint_filter(cache::VariableCache, ::Type{<:MOI.ZeroOne})
    return cache.variable_type == BINARY
end

function check_conflicting_constraint(cache::VariableCache, set::MOI.ZeroOne)
    if !(cache.variable_type == CONTINUOUS)
        throw_conflicting_variable_type(MOI.ZeroOne, cache.variable_type)
    end
    return
end

###
### SingleVariable in Integer
###
### Conflicts with ZeroOne, Semicontinuous, and Semiinteger constraints.

function MOI.add_constraint(model::LinQuadOptimizer, variable::SinVar, set::MOI.Integer)
    __assert_supported_constraint__(model, SinVar, MOI.Integer)
    cache = variable_cache(model, variable)
    check_conflicting_constraint(cache, set)
    cache.variable_type = INTEGER
    change_variable_types!(model, [cache.column], [backend_type(model, set)])
    make_problem_type_integer(model)
    return constraint_index(variable, set)
end

function MOI.delete(model::LinQuadOptimizer, index::SVCI{MOI.Integer})
    __assert_valid__(model, index)
    cache = variable_cache(model, index)
    cache.variable_type = CONTINUOUS
    change_variable_types!(
        model, [cache.column], [backend_type(model, Val{:Continuous}())])
    if !has_integer(model)
        make_problem_type_continuous(model)
    end
end

function MOI.get(::LinQuadOptimizer, ::MOI.ConstraintSet, ::SVCI{MOI.Integer})
    return MOI.Integer()
end

function has_constraint_filter(cache::VariableCache, ::Type{<:MOI.Integer})
    return cache.variable_type == INTEGER
end

function check_conflicting_constraint(cache::VariableCache, set::MOI.Integer)
    if !(cache.variable_type == CONTINUOUS)
        throw_conflicting_variable_type(MOI.Integer, cache.variable_type)
    end
    return
end

###
### SingleVariable in Semicontinuous, Semiinteger
###
### These constraints conflict with all other SingleVariable-in-X constraints.

const SEMI_TYPES = Union{MOI.Semicontinuous{Float64}, MOI.Semiinteger{Float64}}

variable_type(::MOI.Semicontinuous) = SEMI_CONTINUOUS
variable_type(::MOI.Semiinteger) = SEMI_INTEGER

function set_cache(cache::VariableCache, set::SEMI_TYPES)
    cache.lower = set.lower
    cache.upper = set.upper
    cache.bound_type = FREE
    cache.variable_type = variable_type(set)
    return
end

function unset_cache(cache::VariableCache, ::Type{<:SEMI_TYPES})
    cache.lower = -Inf
    cache.upper = Inf
    cache.variable_type = CONTINUOUS
    return
end

function MOI.add_constraint(model::LinQuadOptimizer, variable::SinVar, set::SEMI_TYPES)
    __assert_supported_constraint__(model, SinVar, typeof(set))
    cache = variable_cache(model, variable)
    check_conflicting_constraint(cache, set)
    set_cache(cache, set)
    change_variable_types!(model, [cache.column], [backend_type(model, set)])
    change_variable_bounds!(
        model,
        [cache.column, cache.column],
        [set.lower, set.upper],
        [backend_type(model, Val{:Lowerbound}()), backend_type(model, Val{:Upperbound}())]
    )
    make_problem_type_integer(model)
    return constraint_index(variable, set)
end

function MOI.delete(model::LinQuadOptimizer, index::SVCI{S}) where {S<:SEMI_TYPES}
    __assert_valid__(model, index)
    cache = variable_cache(model, index)
    unset_cache(cache, S)
    change_variable_types!(model, [cache.column], [backend_type(model, Val{:Continuous}())])
    change_variable_bounds!(
        model,
        [cache.column, cache.column],
        [-Inf, Inf],
        [backend_type(model, Val{:Lowerbound}()), backend_type(model, Val{:Upperbound}())]
    )
    if !has_integer(model)
        make_problem_type_continuous(model)
    end
end

function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintSet, index::SVCI{S}) where {S<:SEMI_TYPES}
    cache = variable_cache(model, index)
    return S(cache.lower, cache.upper)
end

function has_constraint_filter(cache::VariableCache, ::Type{<:MOI.Semicontinuous})
    return cache.variable_type == SEMI_CONTINUOUS
end

function has_constraint_filter(cache::VariableCache, ::Type{<:MOI.Semiinteger})
    return cache.variable_type == SEMI_INTEGER
end

function check_conflicting_constraint(cache::VariableCache, set::SEMI_TYPES)
    if !(cache.variable_type == CONTINUOUS)
        throw_conflicting_variable_type(typeof(set), cache.variable_type)
    elseif cache.lower != -Inf
        throw_conflicting_lower_bound(typeof(set))
    elseif cache.upper != Inf
        throw_conflicting_upper_bound(typeof(set))
    end
    return
end
