#=
    Variable bounds
        SingleVariable -in- LessThan
        SingleVariable -in- GreaterThan
        SingleVariable -in- EqualTo
        SingleVariable -in- Interval

        SingleVariable -in- ZeroOne
        SingleVariable -in- Integer
        SingleVariable -in- Semiinteger
        SingleVariable -in- Semicontinuous
=#
constrdict(model::LinQuadOptimizer, ::SVCI{LE}) = cmap(model).upper_bound
constrdict(model::LinQuadOptimizer, ::SVCI{GE}) = cmap(model).lower_bound
constrdict(model::LinQuadOptimizer, ::SVCI{EQ}) = cmap(model).fixed_bound
constrdict(model::LinQuadOptimizer, ::SVCI{IV}) = cmap(model).interval_bound

constrdict(model::LinQuadOptimizer, ::SVCI{MOI.ZeroOne}) = cmap(model).binary
constrdict(model::LinQuadOptimizer, ::SVCI{MOI.Integer}) = cmap(model).integer

constrdict(model::LinQuadOptimizer, ::SVCI{MOI.Semicontinuous{Float64}}) = cmap(model).semicontinuous
constrdict(model::LinQuadOptimizer, ::SVCI{MOI.Semiinteger{Float64}}) = cmap(model).semiinteger

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

"""
    has_value(dict::Dict{K, V}, value::V) where {K, V}

Return true if `dict` has a `key=>value` pair with `value`.
"""
function has_value(dict::Dict{K, V}, value::V) where {K, V}
    return value in values(dict)
end

function MOI.add_constraint(model::LinQuadOptimizer, variable::SinVar, set::S) where S <: LinSets
    __assert_supported_constraint__(model, SinVar, S)
    variable_type = model.variable_type[variable.variable]
    if !(variable_type == Continuous || variable_type == Integer)
        error("Cannot set bounds because variable is of type: $(variable_type).")
    end
    set_variable_bound(model, variable, set)
    model.last_constraint_reference += 1
    index = SVCI{S}(model.last_constraint_reference)
    dict = constrdict(model, index)
    dict[index] = variable.variable
    return index
end

function MOI.add_constraints(model::LinQuadOptimizer, variables::Vector{SinVar},
                             sets::Vector{S}) where S <: LinSets
    __assert_supported_constraint__(model, SinVar, S)
    for variable in variables
        variable_type = model.variable_type[variable.variable]
        if !(variable_type == Continuous || variable_type == Integer)
            error("Cannot set bounds because variable is of type: $(variable_type).")
        end
    end
    set_variable_bounds(model, variables, sets)
    indices = SVCI{S}[]
    for variable in variables
        model.last_constraint_reference += 1
        index = SVCI{S}(model.last_constraint_reference)
        dict = constrdict(model, index)
        dict[index] = variable.variable
        push!(indices, index)
    end
    return indices
end

function MOI.delete(model::LinQuadOptimizer, index::SVCI{S}) where S <: LinSets
    __assert_valid__(model, index)
    delete_constraint_name(model, index)
    dict = constrdict(model, index)
    variable = dict[index]
    model.variable_type[variable] = Continuous
    set_variable_bound(model, SinVar(variable), IV(-Inf, Inf))
    delete!(dict, index)
    return
end

# constraint set
function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintSet, index::SVCI{LE})
    MOI.LessThan{Float64}(
        get_variable_upperbound(model, get_column(model, model[index]))
    )
end
function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintSet, index::SVCI{GE})
    MOI.GreaterThan{Float64}(
        get_variable_lowerbound(model, get_column(model, model[index]))
    )
end
function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintSet, index::SVCI{EQ})
    MOI.EqualTo{Float64}(
        get_variable_lowerbound(model, get_column(model, model[index]))
    )
end
function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintSet, index::SVCI{IV})
    column = get_column(model, model[index])
    return MOI.Interval{Float64}(get_variable_lowerbound(model, column),
                                 get_variable_upperbound(model, column))
end

# constraint function
function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintFunction, index::SVCI{<: LinSets})
    return SinVar(model[index])
end

# modify
MOI.supports(::LinQuadOptimizer, ::MOI.ConstraintSet, ::Type{SVCI{S}}) where S<:LinSets = true
function MOI.set(model::LinQuadOptimizer, ::MOI.ConstraintSet, index::SVCI{S}, newset::S) where S <: LinSets
    set_variable_bound(model, SinVar(model[index]), newset)
end

#=
    Binary constraints

For some reason CPLEX doesn't respect bounds on a binary variable, so we should
store the previous bounds so that if we delete the binary constraint we can
revert to the old bounds.

Xpress is worse, once binary, the bounds are changed independently of what the
user does.
=#
function MOI.add_constraint(model::LinQuadOptimizer, variable::SinVar, set::MOI.ZeroOne)
    __assert_supported_constraint__(model, SinVar, MOI.ZeroOne)
    variable_type = model.variable_type[variable.variable]
    if variable_type != Continuous
        error("Cannot make variable binary because it is $(variable_type).")
    end
    model.variable_type[variable.variable] = Binary
    model.last_constraint_reference += 1
    index = SVCI{MOI.ZeroOne}(model.last_constraint_reference)
    column = get_column(model, variable)
    dict = constrdict(model, index)
    dict[index] = (variable.variable, get_variable_lowerbound(model, column),
                   get_variable_upperbound(model, column))
    change_variable_types!(model, [column], [backend_type(model, set)])
    change_variable_bounds!(model,
        [column, column],
        [0.0, 1.0],
        [backend_type(model, Val{:Lowerbound}()),
         backend_type(model, Val{:Upperbound}())]
    )
    make_problem_type_integer(model)
    return index
end

function MOI.delete(model::LinQuadOptimizer, index::SVCI{MOI.ZeroOne})
    __assert_valid__(model, index)
    delete_constraint_name(model, index)
    dict = constrdict(model, index)
    (variable, lower, upper) = dict[index]
    model.variable_type[variable] = Continuous
    column = get_column(model, variable)
    change_variable_types!(
        model, [column], [backend_type(model, Val{:Continuous}())])
    change_variable_bounds!(model,
       [column, column],
       [lower, upper],
       [backend_type(model, Val{:Lowerbound}()),
        backend_type(model, Val{:Upperbound}())]
    )
    delete!(dict, index)
    if !has_integer(model)
        make_problem_type_continuous(model)
    end
end

function MOI.get(::LinQuadOptimizer, ::MOI.ConstraintSet, ::SVCI{MOI.ZeroOne})
    return MOI.ZeroOne()
end

function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintFunction, index::SVCI{MOI.ZeroOne})
    return SinVar(model[index][1])
end

#=
    Integer constraints
=#

function MOI.add_constraint(model::LinQuadOptimizer, variable::SinVar, set::MOI.Integer)
    __assert_supported_constraint__(model, SinVar, MOI.Integer)
    variable_type = model.variable_type[variable.variable]
    if variable_type != Continuous
        error("Cannot make variable integer because it is $(variable_type).")
    end
    model.variable_type[variable.variable] = Integer
    change_variable_types!(model, [get_column(model, variable)],
                           [backend_type(model, set)])
    model.last_constraint_reference += 1
    index = SVCI{MOI.Integer}(model.last_constraint_reference)
    dict = constrdict(model, index)
    dict[index] = variable.variable
    make_problem_type_integer(model)
    return index
end

function MOI.delete(model::LinQuadOptimizer, index::SVCI{MOI.Integer})
    __assert_valid__(model, index)
    delete_constraint_name(model, index)
    dict = constrdict(model, index)
    variable = dict[index]
    model.variable_type[variable] = Continuous
    change_variable_types!(model, [get_column(model, variable)],
                           [backend_type(model, Val{:Continuous}())])
    delete!(dict, index)
    if !has_integer(model)
        make_problem_type_continuous(model)
    end
end

function MOI.get(::LinQuadOptimizer, ::MOI.ConstraintSet, ::SVCI{MOI.Integer})
    return MOI.Integer()
end

function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintFunction,
                 index::SVCI{MOI.Integer})
    return SinVar(model[index])
end

#=
    Semicontinuous / Semiinteger constraints
=#
const SEMI_TYPES = Union{MOI.Semicontinuous{Float64}, MOI.Semiinteger{Float64}}
variable_type_(::Type{<:MOI.Semicontinuous}) = Semicontinuous
variable_type_(::Type{<:MOI.Semiinteger}) = Semiinteger
function MOI.add_constraint(model::LinQuadOptimizer, variable::SinVar, set::S) where S <: SEMI_TYPES
    __assert_supported_constraint__(model, SinVar, S)
    variable_type = model.variable_type[variable.variable]
    if variable_type != Continuous
        error("Cannot make variable $(S) because it is $(variable_type).")
    end
    model.variable_type[variable.variable] = variable_type_(S)
    column = get_column(model, variable)
    change_variable_types!(model, [column], [backend_type(model, set)])
    change_variable_bounds!(model,
        [column, column],
        [set.lower, set.upper],
        [backend_type(model, Val{:Lowerbound}()),
         backend_type(model, Val{:Upperbound}())]
    )
    model.last_constraint_reference += 1
    index = SVCI{S}(model.last_constraint_reference)
    dict = constrdict(model, index)
    dict[index] = variable.variable
    make_problem_type_integer(model)
    return index
end

function MOI.delete(model::LinQuadOptimizer, index::SVCI{<:SEMI_TYPES})
    __assert_valid__(model, index)
    delete_constraint_name(model, index)
    dict = constrdict(model, index)
    variable = dict[index]
    column = get_column(model, variable)
    model.variable_type[variable] = Continuous
    change_variable_types!(model, [column], [backend_type(model, Val{:Continuous}())])
    change_variable_bounds!(model,
        [column, column],
        [-Inf, Inf],
        [backend_type(model, Val{:Lowerbound}()),
         backend_type(model, Val{:Upperbound}())]
    )
    delete!(dict, index)
    if !has_integer(model)
        make_problem_type_continuous(model)
    end
end

function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintSet, index::SVCI{S}) where S <: SEMI_TYPES
    dict = constrdict(model, index)
    column = get_column(model, dict[index])
    return S(get_variable_lowerbound(model, column),
             get_variable_upperbound(model, column))
end

function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintFunction, index::SVCI{<:SEMI_TYPES})
    return SinVar(model[index])
end
