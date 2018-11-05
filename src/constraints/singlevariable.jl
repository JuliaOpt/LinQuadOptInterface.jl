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

function set_variable_bound(model::LinQuadOptimizer, v::SinVar, set::LE)
    change_variable_bounds!(model,
        [get_column(model, v)],
        [set.upper],
        [backend_type(model, Val{:Upperbound}())]
    )
end

function set_variable_bound(model::LinQuadOptimizer, v::SinVar, set::GE)
    change_variable_bounds!(model,
        [get_column(model, v)],
        [set.lower],
        [backend_type(model, Val{:Lowerbound}())]
    )
end

function set_variable_bound(model::LinQuadOptimizer, v::SinVar, set::EQ)
    change_variable_bounds!(model,
        [get_column(model, v), get_column(model, v)],
        [set.value, set.value],
        [backend_type(model, Val{:Upperbound}()),
         backend_type(model, Val{:Lowerbound}())
        ]
    )
end

function set_variable_bound(model::LinQuadOptimizer, v::SinVar, set::IV)
    change_variable_bounds!(model,
        [get_column(model, v), get_column(model, v)],
        [set.upper, set.lower],
        [backend_type(model, Val{:Upperbound}()),
         backend_type(model, Val{:Lowerbound}())
        ]
    )
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
