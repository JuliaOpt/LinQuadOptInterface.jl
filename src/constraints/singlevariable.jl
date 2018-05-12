#=
    Variable bounds
        SingleVariable -in- LessThan
        SingleVariable -in- GreaterThan
        SingleVariable -in- EqualTo
        SingleVariable -in- Interval

TODO

    Binary
        SingleVariable -in- ZeroOne
        SingleVariable -in- Integer
=#
constrdict(m::LinQuadOptimizer, ::SVCI{LE}) = cmap(m).upper_bound
constrdict(m::LinQuadOptimizer, ::SVCI{GE}) = cmap(m).lower_bound
constrdict(m::LinQuadOptimizer, ::SVCI{EQ}) = cmap(m).fixed_bound
constrdict(m::LinQuadOptimizer, ::SVCI{IV}) = cmap(m).interval_bound

constrdict(m::LinQuadOptimizer, ::SVCI{MOI.ZeroOne}) = cmap(m).binary
constrdict(m::LinQuadOptimizer, ::SVCI{MOI.Integer}) = cmap(m).integer

function setvariablebound!(m::LinQuadOptimizer, col::Int, bound::Float64, sense::Cchar)
    change_variable_bounds!(m, [col], [bound], [sense])
end

function setvariablebound!(m::LinQuadOptimizer, v::SinVar, set::LE)
    setvariablebound!(m, getcol(m, v), set.upper, backend_type(m, Val{:Upperbound}()))
end
function setvariablebound!(m::LinQuadOptimizer, v::SinVar, set::GE)
    setvariablebound!(m, getcol(m, v), set.lower, backend_type(m, Val{:Lowerbound}()))
end
function setvariablebound!(m::LinQuadOptimizer, v::SinVar, set::EQ)
    setvariablebound!(m, getcol(m, v), set.value, backend_type(m, Val{:Upperbound}()))
    setvariablebound!(m, getcol(m, v), set.value, backend_type(m, Val{:Lowerbound}()))
end
function setvariablebound!(m::LinQuadOptimizer, v::SinVar, set::IV)
    setvariablebound!(m, getcol(m, v), set.upper, backend_type(m, Val{:Upperbound}()))
    setvariablebound!(m, getcol(m, v), set.lower, backend_type(m, Val{:Lowerbound}()))
end

# add constraint
function MOI.addconstraint!(m::LinQuadOptimizer, v::SinVar, set::S) where S <: LinSets
    setvariablebound!(m, v, set)
    m.last_constraint_reference += 1
    ref = SVCI{S}(m.last_constraint_reference)
    dict = constrdict(m, ref)
    dict[ref] = v.variable
    ref
end

# delete constraint
MOI.candelete(m::LinQuadOptimizer, c::SVCI{S}) where S <: LinSets = MOI.isvalid(m, c)
function MOI.delete!(m::LinQuadOptimizer, c::SVCI{S}) where S <: LinSets
    deleteconstraintname!(m, c)
    dict = constrdict(m, c)
    vref = dict[c]
    setvariablebound!(m, SinVar(vref), IV(-Inf, Inf))
    delete!(dict, c)
end

# constraint set
MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintSet, ::Type{SVCI{S}}) where S <: LinSets = true
function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintSet, c::SVCI{LE})
    MOI.LessThan{Float64}(
        get_variable_upperbound(m, getcol(m, m[c]))
    )
end
function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintSet, c::SVCI{GE})
    MOI.GreaterThan{Float64}(
        get_variable_lowerbound(m, getcol(m, m[c]))
    )
end
function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintSet, c::SVCI{EQ})
    MOI.EqualTo{Float64}(
        get_variable_lowerbound(m, getcol(m, m[c]))
    )
end
function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintSet, c::SVCI{IV})
    lb = get_variable_lowerbound(m, getcol(m, m[c]))
    ub = get_variable_upperbound(m, getcol(m, m[c]))
    return MOI.Interval{Float64}(lb, ub)
end

# constraint function
MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintFunction, ::Type{SVCI{S}}) where S <: LinSets = true
function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintFunction, c::SVCI{<: LinSets})
    return SinVar(m[c])
end

# modify
MOI.canmodifyconstraint(::LinQuadOptimizer, ::SVCI{S}, ::Type{S}) where S <: LinSets = true
function MOI.modifyconstraint!(m::LinQuadOptimizer, c::SVCI{S}, newset::S) where S <: LinSets
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
    m.last_constraint_reference += 1
    ref = SVCI{MOI.ZeroOne}(m.last_constraint_reference)
    dict = constrdict(m, ref)
    ub = get_variable_upperbound(m, getcol(m, v))
    lb = get_variable_lowerbound(m, getcol(m, v))
    dict[ref] = (v.variable, lb, ub)
    change_variable_types!(m, [getcol(m, v)], [backend_type(m, set)])
    setvariablebound!(m, getcol(m, v), 1.0, backend_type(m, Val{:Upperbound}()))
    setvariablebound!(m, getcol(m, v), 0.0, backend_type(m, Val{:Lowerbound}()))
    make_problem_type_integer(m)
    ref
end

MOI.candelete(m::LinQuadOptimizer, c::SVCI{MOI.ZeroOne}) = MOI.isvalid(m, c)
function MOI.delete!(m::LinQuadOptimizer, c::SVCI{MOI.ZeroOne})
    deleteconstraintname!(m, c)
    dict = constrdict(m, c)
    (v, lb, ub) = dict[c]
    change_variable_types!(m, [getcol(m, v)], [backend_type(m, Val{:Continuous}())])
    setvariablebound!(m, getcol(m, v), ub, backend_type(m, Val{:Upperbound}()))
    setvariablebound!(m, getcol(m, v), lb, backend_type(m, Val{:Lowerbound}()))
    delete!(dict, c)
    if !hasinteger(m)
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
    change_variable_types!(m, [getcol(m, v)], [backend_type(m, set)])
    m.last_constraint_reference += 1
    ref = SVCI{MOI.Integer}(m.last_constraint_reference)
    dict = constrdict(m, ref)
    dict[ref] = v.variable
    make_problem_type_integer(m)
    ref
end

MOI.candelete(m::LinQuadOptimizer, c::SVCI{MOI.Integer}) = MOI.isvalid(m, c)
function MOI.delete!(m::LinQuadOptimizer, c::SVCI{MOI.Integer})
    deleteconstraintname!(m, c)
    dict = constrdict(m, c)
    v = dict[c]
    change_variable_types!(m, [getcol(m, v)], [backend_type(m, Val{:Continuous}())])
    delete!(dict, c)
    if !hasinteger(m)
        make_problem_type_continuous(m)
    end
end

MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintSet, ::Type{SVCI{MOI.Integer}}) = true
MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintSet, c::SVCI{MOI.Integer}) = MOI.Integer()

MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintFunction, ::Type{SVCI{MOI.Integer}}) = true
MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintFunction, c::SVCI{MOI.Integer}) = SinVar(m[c])
