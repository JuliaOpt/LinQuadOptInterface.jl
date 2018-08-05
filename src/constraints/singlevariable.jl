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
constrdict(m::LinQuadOptimizer, ::SVCI{LE}) = cmap(m).upper_bound
constrdict(m::LinQuadOptimizer, ::SVCI{GE}) = cmap(m).lower_bound
constrdict(m::LinQuadOptimizer, ::SVCI{EQ}) = cmap(m).fixed_bound
constrdict(m::LinQuadOptimizer, ::SVCI{IV}) = cmap(m).interval_bound

constrdict(m::LinQuadOptimizer, ::SVCI{MOI.ZeroOne}) = cmap(m).binary
constrdict(m::LinQuadOptimizer, ::SVCI{MOI.Integer}) = cmap(m).integer

constrdict(m::LinQuadOptimizer, ::SVCI{MOI.Semicontinuous{Float64}}) = cmap(m).semicontinuous
constrdict(m::LinQuadOptimizer, ::SVCI{MOI.Semiinteger{Float64}}) = cmap(m).semiinteger

function setvariablebound!(m::LinQuadOptimizer, col::Int, bound::Float64, sense::Cchar)
    change_variable_bounds!(m, [col], [bound], [sense])
end

function setvariablebound!(m::LinQuadOptimizer, v::SinVar, set::LE)
    setvariablebound!(m, get_column(m, v), set.upper, backend_type(m, Val{:Upperbound}()))
end
function setvariablebound!(m::LinQuadOptimizer, v::SinVar, set::GE)
    setvariablebound!(m, get_column(m, v), set.lower, backend_type(m, Val{:Lowerbound}()))
end
function setvariablebound!(m::LinQuadOptimizer, v::SinVar, set::EQ)
    setvariablebound!(m, get_column(m, v), set.value, backend_type(m, Val{:Upperbound}()))
    setvariablebound!(m, get_column(m, v), set.value, backend_type(m, Val{:Lowerbound}()))
end
function setvariablebound!(m::LinQuadOptimizer, v::SinVar, set::IV)
    setvariablebound!(m, get_column(m, v), set.upper, backend_type(m, Val{:Upperbound}()))
    setvariablebound!(m, get_column(m, v), set.lower, backend_type(m, Val{:Lowerbound}()))
end

SVCI(v::SinVar, ::S) where S = SVCI{S}(v.variable.value)

function hasvalue(d::Dict, val)
    for v in values(d)
        if v == val
            return true
        end
    end
    return false
end

function checkexisting(m::LinQuadOptimizer, v::SinVar, set::S) where S
    ref = SVCI(v, set)
    if hasvalue(constrdict(m, ref), v.variable)
        error("Adding the same constraint type: $(S) is not allowed for SingleVariable function")
    end
end

function checkconflicting(m::LinQuadOptimizer, v::SinVar, set_to_add::S0, set_to_test::S) where S where S0
    ref = SVCI(v, set_to_test)
    if hasvalue(constrdict(m, ref), v.variable)
        error("Adding the same constraint type: $(S0) is not allowed for SingleVariable function because there is constraint of type $(S) tied to the respective variable")
    end
end

# add constraint
function MOI.addconstraint!(m::LinQuadOptimizer, v::SinVar, set::S) where S <: LinSets
    __assert_supported_constraint__(m, SinVar, S)
    checkexisting(m, v, set)
    checkconflicting(m, v, set, MOI.Semicontinuous(0.0, 0.0))
    checkconflicting(m, v, set, MOI.Semiinteger(0.0, 0.0))
    checkconflicting(m, v, set, MOI.ZeroOne())
    setvariablebound!(m, v, set)
    m.last_constraint_reference += 1
    ref = SVCI{S}(m.last_constraint_reference)
    # ref = SVCI(v, set)
    dict = constrdict(m, ref)
    dict[ref] = v.variable
    ref
end

# delete constraint
function MOI.delete!(m::LinQuadOptimizer, c::SVCI{S}) where S <: LinSets
    __assert_valid__(m, c)
    delete_constraint_name(m, c)
    dict = constrdict(m, c)
    vref = dict[c]
    setvariablebound!(m, SinVar(vref), IV(-Inf, Inf))
    delete!(dict, c)
end

# constraint set
MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintSet, ::Type{SVCI{S}}) where S <: LinSets = true
function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintSet, c::SVCI{LE})
    MOI.LessThan{Float64}(
        get_variable_upperbound(m, get_column(m, m[c]))
    )
end
function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintSet, c::SVCI{GE})
    MOI.GreaterThan{Float64}(
        get_variable_lowerbound(m, get_column(m, m[c]))
    )
end
function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintSet, c::SVCI{EQ})
    MOI.EqualTo{Float64}(
        get_variable_lowerbound(m, get_column(m, m[c]))
    )
end
function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintSet, c::SVCI{IV})
    lb = get_variable_lowerbound(m, get_column(m, m[c]))
    ub = get_variable_upperbound(m, get_column(m, m[c]))
    return MOI.Interval{Float64}(lb, ub)
end

# constraint function
MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintFunction, ::Type{SVCI{S}}) where S <: LinSets = true
function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintFunction, c::SVCI{<: LinSets})
    return SinVar(m[c])
end

# modify
MOI.supports(::LinQuadOptimizer, ::MOI.ConstraintSet, ::Type{SVCI{S}}) where S<:LinSets = true
function MOI.set!(m::LinQuadOptimizer, ::MOI.ConstraintSet, c::SVCI{S}, newset::S) where S <: LinSets
    setvariablebound!(m, SinVar(m[c]), newset)
end

#=
    Binary constraints

For some reason CPLEX doesn't respect bounds on a binary variable, so we
should store the previous bounds so that if we delete the binary constraint
we can revert to the old bounds

Xpress is worse, once binary, the bounds are changed independently of what the user does
=#
function MOI.addconstraint!(m::LinQuadOptimizer, v::SinVar, set::MOI.ZeroOne)
    __assert_supported_constraint__(m, SinVar, MOI.ZeroOne)
    checkexisting(m, v, set)
    checkconflicting(m, v, set, MOI.Integer())
    checkconflicting(m, v, set, MOI.Semicontinuous(0.0, 0.0))
    checkconflicting(m, v, set, MOI.Semiinteger(0.0, 0.0))
    m.last_constraint_reference += 1
    ref = SVCI{MOI.ZeroOne}(m.last_constraint_reference)
    # ref = SVCI(v, set)
    dict = constrdict(m, ref)
    ub = get_variable_upperbound(m, get_column(m, v))
    lb = get_variable_lowerbound(m, get_column(m, v))
    dict[ref] = (v.variable, lb, ub)
    change_variable_types!(m, [get_column(m, v)], [backend_type(m, set)])
    setvariablebound!(m, get_column(m, v), 1.0, backend_type(m, Val{:Upperbound}()))
    setvariablebound!(m, get_column(m, v), 0.0, backend_type(m, Val{:Lowerbound}()))
    make_problem_type_integer(m)
    ref
end

function MOI.delete!(m::LinQuadOptimizer, c::SVCI{MOI.ZeroOne})
    __assert_valid__(m, c)
    delete_constraint_name(m, c)
    dict = constrdict(m, c)
    (v, lb, ub) = dict[c]
    change_variable_types!(m, [get_column(m, v)], [backend_type(m, Val{:Continuous}())])
    setvariablebound!(m, get_column(m, v), ub, backend_type(m, Val{:Upperbound}()))
    setvariablebound!(m, get_column(m, v), lb, backend_type(m, Val{:Lowerbound}()))
    delete!(dict, c)
    if !has_integer(m)
        make_problem_type_continuous(m)
    end
end

MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintSet, ::Type{SVCI{MOI.ZeroOne}}) = true
MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintSet, c::SVCI{MOI.ZeroOne}) = MOI.ZeroOne()

MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintFunction, ::Type{SVCI{MOI.ZeroOne}}) = true
MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintFunction, c::SVCI{MOI.ZeroOne}) = SinVar(m[c][1])

#=
    Integer constraints
=#

function MOI.addconstraint!(m::LinQuadOptimizer, v::SinVar, set::MOI.Integer)
    __assert_supported_constraint__(m, SinVar, MOI.Integer)
    checkexisting(m, v, set)
    checkconflicting(m, v, set, MOI.ZeroOne())
    checkconflicting(m, v, set, MOI.Semicontinuous(0.0, 0.0))
    checkconflicting(m, v, set, MOI.Semiinteger(0.0, 0.0))
    change_variable_types!(m, [get_column(m, v)], [backend_type(m, set)])
    m.last_constraint_reference += 1
    ref = SVCI{MOI.Integer}(m.last_constraint_reference)
    # ref = SVCI(v, set)
    dict = constrdict(m, ref)
    dict[ref] = v.variable
    make_problem_type_integer(m)
    ref
end

function MOI.delete!(m::LinQuadOptimizer, c::SVCI{MOI.Integer})
    __assert_valid__(m, c)
    delete_constraint_name(m, c)
    dict = constrdict(m, c)
    v = dict[c]
    change_variable_types!(m, [get_column(m, v)], [backend_type(m, Val{:Continuous}())])
    delete!(dict, c)
    if !has_integer(m)
        make_problem_type_continuous(m)
    end
end

MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintSet, ::Type{SVCI{MOI.Integer}}) = true
MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintSet, c::SVCI{MOI.Integer}) = MOI.Integer()

MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintFunction, ::Type{SVCI{MOI.Integer}}) = true
MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintFunction, c::SVCI{MOI.Integer}) = SinVar(m[c])

#=
    Semicontinuous / Semiinteger constraints
=#
const SEMI_TYPES = Union{MOI.Semicontinuous{Float64}, MOI.Semiinteger{Float64}}
function MOI.addconstraint!(m::LinQuadOptimizer, v::SinVar, set::S) where S <: SEMI_TYPES
    __assert_supported_constraint__(m, SinVar, S)
    checkexisting(m, v, set)
    checkconflicting(m, v, set, MOI.ZeroOne())
    checkconflicting(m, v, set, MOI.Integer())
    if S == MOI.Semicontinuous{Float64}
        checkconflicting(m, v, set, MOI.Semiinteger(0.0, 0.0))
    else
        checkconflicting(m, v, set, MOI.Semicontinuous(0.0, 0.0))
    end
    change_variable_types!(m, [get_column(m, v)], [backend_type(m, set)])
    setvariablebound!(m, get_column(m, v), set.upper, backend_type(m, Val{:Upperbound}()))
    setvariablebound!(m, get_column(m, v), set.lower, backend_type(m, Val{:Lowerbound}()))
    m.last_constraint_reference += 1
    ref = SVCI{S}(m.last_constraint_reference)
    # ref = SVCI(v, set)
    dict = constrdict(m, ref)
    dict[ref] = v.variable
    make_problem_type_integer(m)
    ref
end

function MOI.delete!(m::LinQuadOptimizer, c::SVCI{<:SEMI_TYPES})
    __assert_valid__(m, c)
    delete_constraint_name(m, c)
    dict = constrdict(m, c)
    v = dict[c]
    change_variable_types!(m, [get_column(m, v)], [backend_type(m, Val{:Continuous}())])
    setvariablebound!(m, get_column(m, v), Inf, backend_type(m, Val{:Upperbound}()))
    setvariablebound!(m, get_column(m, v), -Inf, backend_type(m, Val{:Lowerbound}()))
    delete!(dict, c)
    if !has_integer(m)
        make_problem_type_continuous(m)
    end
end

MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintSet, ::Type{SVCI{S}}) where S <:SEMI_TYPES = true
function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintSet, c::SVCI{S}) where S <: SEMI_TYPES
    dict = constrdict(m, c)
    v = dict[c]
    lb = get_variable_lowerbound(m, get_column(m, v))
    ub = get_variable_upperbound(m, get_column(m, v))
    return S(lb, ub)
end

MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintFunction, ::Type{SVCI{S}}) where S <:SEMI_TYPES = true
MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintFunction, c::SVCI{<:SEMI_TYPES}) = SinVar(m[c])
