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

function set_variable_bound(model::LinQuadOptimizer, column::Int,
                            bound::Float64, sense)
    change_variable_bounds!(model, [column], [bound], [sense])
end

function set_variable_bound(model::LinQuadOptimizer, v::SinVar, set::LE)
    set_variable_bound(model, get_column(model, v), set.upper,
                      backend_type(model, Val{:Upperbound}()))
end
function set_variable_bound(model::LinQuadOptimizer, v::SinVar, set::GE)
    set_variable_bound(model, get_column(model, v), set.lower,
                      backend_type(model, Val{:Lowerbound}()))
end
function set_variable_bound(model::LinQuadOptimizer, v::SinVar, set::EQ)
    set_variable_bound(model, get_column(model, v), set.value,
                      backend_type(model, Val{:Upperbound}()))
    set_variable_bound(model, get_column(model, v), set.value,
                      backend_type(model, Val{:Lowerbound}()))
end
function set_variable_bound(model::LinQuadOptimizer, v::SinVar, set::IV)
    set_variable_bound(model, get_column(model, v), set.upper,
                      backend_type(model, Val{:Upperbound}()))
    set_variable_bound(model, get_column(model, v), set.lower,
                      backend_type(model, Val{:Lowerbound}()))
end

"""
    has_value(dict::Dict{K, V}, value::V) where {K, V}

Return true if `dict` has a `key=>value` pair with `value`.
"""
function has_value(dict::Dict{K, V}, value::V) where {K, V}
    return value in values(dict)
end

"""
    __check_for_conflicting__(model::LinQuadOptimizer, variable::SinVar, set,
                              conflicting_type::Type{<:AbstractSet})

Throw an error if `variable` is already constrained to be in a set of type
`conflicting_type`.
"""
function __check_for_conflicting__(model::LinQuadOptimizer, variable::SinVar, set,
                                   conflict_type::Type{<:MOI.AbstractSet})
    if has_value(constrdict(model, SVCI{conflict_type}(0)), variable.variable)
        error("Cannot add constraint $(variable)-in-$(set) as it is already " *
              "constrained by a set of type $(conflict_type).")
    end
end

function __check_for_conflicting__(model::LinQuadOptimizer, variable::SinVar, set,
                                   conflict_type::Type{MOI.ZeroOne})
    for (index, lower, upper) in values(constrdict(model, SVCI{conflict_type}(0)))
        if index == variable.variable
            error("Cannot add constraint $(variable)-in-$(set) as it is already " *
                  "constrained by a set of type $(conflict_type).")
        end
    end
end


"""
    __check_for_conflicting__(model::LinQuadOptimizer, variable::SinVar, set,
                              conflicting_types::AbstractSet...)

Throw an error if `variable` is constrained to be in a set whose type is one of
`conflicting_types`.
"""
function __check_for_conflicting__(model::LinQuadOptimizer, variable::SinVar, set,
                                   conflicting_types...)
    for conflict_type in conflicting_types
        __check_for_conflicting__(model, variable, set, conflict_type)
    end
end

function MOI.add_constraint(model::LinQuadOptimizer, variable::SinVar, set::S) where S <: LinSets
    __assert_supported_constraint__(model, SinVar, S)
    # Since the following "variable type" sets also define bounds (implicitly or explicitly), 
    # they may conflict with other bound constraints.
    __check_for_conflicting__(model, variable, set,
        S, MOI.Semicontinuous{Float64}, MOI.Semiinteger{Float64}, MOI.ZeroOne)
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
    variable_index = dict[index]
    set_variable_bound(model, SinVar(variable_index), IV(-Inf, Inf))
    delete!(dict, index)
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

For some reason CPLEX doesn't respect bounds on a binary variable, so we
should store the previous bounds so that if we delete the binary constraint
we can revert to the old bounds

Xpress is worse, once binary, the bounds are changed independently of what the user does
=#
function MOI.add_constraint(model::LinQuadOptimizer, variable::SinVar, set::MOI.ZeroOne)
    __assert_supported_constraint__(model, SinVar, MOI.ZeroOne)
    __check_for_conflicting__(model, variable, set, MOI.ZeroOne, MOI.Integer,
        MOI.Semicontinuous{Float64}, MOI.Semiinteger{Float64})
    model.last_constraint_reference += 1
    index = SVCI{MOI.ZeroOne}(model.last_constraint_reference)
    column = get_column(model, variable)
    dict = constrdict(model, index)
    dict[index] = (variable.variable, get_variable_lowerbound(model, column),
                    get_variable_upperbound(model, column))
    change_variable_types!(model, [column], [backend_type(model, set)])
    set_variable_bound(model, column, 1.0,
                       backend_type(model, Val{:Upperbound}()))
    set_variable_bound(model, column, 0.0,
                       backend_type(model, Val{:Lowerbound}()))
    make_problem_type_integer(model)
    return index
end

function MOI.delete(model::LinQuadOptimizer, index::SVCI{MOI.ZeroOne})
    __assert_valid__(model, index)
    delete_constraint_name(model, index)
    dict = constrdict(model, index)
    (variable, lower, upper) = dict[index]
    column = get_column(model, variable)
    change_variable_types!(model, [column],
                           [backend_type(model, Val{:Continuous}())])
    set_variable_bound(model, column, upper,
                       backend_type(model, Val{:Upperbound}()))
    set_variable_bound(model, column, lower,
                       backend_type(model, Val{:Lowerbound}()))
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
    __check_for_conflicting__(model, variable, set, MOI.ZeroOne,
        MOI.Semicontinuous{Float64}, MOI.Semiinteger{Float64})
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
function MOI.add_constraint(model::LinQuadOptimizer, variable::SinVar, set::S) where S <: SEMI_TYPES
    __assert_supported_constraint__(model, SinVar, S)
    __check_for_conflicting__(model, variable, set, S, MOI.ZeroOne, MOI.Integer)
    if S == MOI.Semicontinuous{Float64}
        __check_for_conflicting__(model, variable, set, MOI.Semiinteger{Float64})
    else
        __check_for_conflicting__(model, variable, set, MOI.Semicontinuous{Float64})
    end
    column = get_column(model, variable)
    change_variable_types!(model, [column], [backend_type(model, set)])
    set_variable_bound(model, column, set.upper,
                       backend_type(model, Val{:Upperbound}()))
    set_variable_bound(model, column, set.lower,
                       backend_type(model, Val{:Lowerbound}()))
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
    column = get_column(model, dict[index])
    change_variable_types!(model, [column], [backend_type(model, Val{:Continuous}())])
    set_variable_bound(model, column, Inf,
                       backend_type(model, Val{:Upperbound}()))
    set_variable_bound(model, column, -Inf,
                       backend_type(model, Val{:Lowerbound}()))
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
