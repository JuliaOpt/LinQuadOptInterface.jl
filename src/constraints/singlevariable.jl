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
    lqs_chgbds!(m, [col], [bound], [sense])
end

function setvariablebound!(m::LinQuadOptimizer, v::SinVar, set::LE)
    setvariablebound!(m, getcol(m, v), set.upper, lqs_char(m, Val{:Upperbound}()))
end
function setvariablebound!(m::LinQuadOptimizer, v::SinVar, set::GE)
    setvariablebound!(m, getcol(m, v), set.lower, lqs_char(m, Val{:Lowerbound}()))
end
function setvariablebound!(m::LinQuadOptimizer, v::SinVar, set::EQ)
    setvariablebound!(m, getcol(m, v), set.value, lqs_char(m, Val{:Upperbound}()))
    setvariablebound!(m, getcol(m, v), set.value, lqs_char(m, Val{:Lowerbound}()))
end
function setvariablebound!(m::LinQuadOptimizer, v::SinVar, set::IV)
    setvariablebound!(m, getcol(m, v), set.upper, lqs_char(m, Val{:Upperbound}()))
    setvariablebound!(m, getcol(m, v), set.lower, lqs_char(m, Val{:Lowerbound}()))
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
MOI.candelete(m::LinQuadOptimizer, c::SVCI{S}) where S <: LinSets = true
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
    MOI.LessThan{Float64}(lqs_getub(m, getcol(m, m[c])))
end
function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintSet, c::SVCI{GE})
    MOI.GreaterThan{Float64}(lqs_getlb(m, getcol(m, m[c])))
end
function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintSet, c::SVCI{EQ})
    MOI.EqualTo{Float64}(lqs_getlb(m, getcol(m, m[c])))
end
function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintSet, c::SVCI{IV})
    lb = lqs_getlb(m, getcol(m, m[c]))
    ub = lqs_getub(m, getcol(m, m[c]))
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
    ub = lqs_getub(m, getcol(m, v))
    lb = lqs_getlb(m, getcol(m, v))
    dict[ref] = (v.variable, lb, ub)
    lqs_chgctype!(m, [getcol(m, v)], [lqs_char(m, set)])
    setvariablebound!(m, getcol(m, v), 1.0, lqs_char(m, Val{:Upperbound}()))
    setvariablebound!(m, getcol(m, v), 0.0, lqs_char(m, Val{:Lowerbound}()))
    lqs_make_problem_type_integer(m)
    ref
end

MOI.candelete(m::LinQuadOptimizer, c::SVCI{MOI.ZeroOne}) = true
function MOI.delete!(m::LinQuadOptimizer, c::SVCI{MOI.ZeroOne})
    deleteconstraintname!(m, c)
    dict = constrdict(m, c)
    (v, lb, ub) = dict[c]
    lqs_chgctype!(m, [getcol(m, v)], [lqs_char(m, Val{:Continuous}())])
    setvariablebound!(m, getcol(m, v), ub, lqs_char(m, Val{:Upperbound}()))
    setvariablebound!(m, getcol(m, v), lb, lqs_char(m, Val{:Lowerbound}()))
    delete!(dict, c)
    if !hasinteger(m)
        lqs_make_problem_type_continuous(m)
    end
end

MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintSet, ::Type{SVCI{MOI.ZeroOne}}) = true
MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintSet, c::SVCI{MOI.ZeroOne}) = MOI.ZeroOne()

MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintFunction, ::Type{SVCI{MOI.ZeroOne}}) = true
MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintFunction, c::SVCI{MOI.ZeroOne}) = m[c]

#=
    Integer constraints
=#

function MOI.addconstraint!(m::LinQuadOptimizer, v::SinVar, set::MOI.Integer)
    lqs_chgctype!(m, [getcol(m, v)], [lqs_char(m, set)])
    m.last_constraint_reference += 1
    ref = SVCI{MOI.Integer}(m.last_constraint_reference)
    dict = constrdict(m, ref)
    dict[ref] = v.variable
    lqs_make_problem_type_integer(m)
    ref
end

function MOI.delete!(m::LinQuadOptimizer, c::SVCI{MOI.Integer})
    deleteconstraintname!(m, c)
    dict = constrdict(m, c)
    v = dict[c]
    lqs_chgctype!(m, [getcol(m, v)], [lqs_char(m, Val{:Continuous}())])
    delete!(dict, c)
    if !hasinteger(m)
        lqs_make_problem_type_continuous(m)
    end
end
MOI.candelete(m::LinQuadOptimizer, c::SVCI{MOI.Integer}) = true

MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintSet, ::Type{SVCI{MOI.Integer}}) = true
MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintSet, c::SVCI{MOI.Integer}) = MOI.Integer()

MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintFunction, c::SVCI{MOI.Integer}) = m[c]
MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintFunction, ::Type{SVCI{MOI.Integer}}) = true
