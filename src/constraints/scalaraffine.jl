#=
        ScalarAffineFunction -in- LessThan
        ScalarAffineFunction -in- GreaterThan
        ScalarAffineFunction -in- EqualTo
        ScalarAffineFunction -in- Interval
=#
constrdict(m::LinQuadOptimizer, ::LCI{LE})  = cmap(m).less_than
constrdict(m::LinQuadOptimizer, ::LCI{GE})  = cmap(m).greater_than
constrdict(m::LinQuadOptimizer, ::LCI{EQ})  = cmap(m).equal_to
constrdict(m::LinQuadOptimizer, ::LCI{IV})  = cmap(m).interval

function MOI.addconstraint!(m::LinQuadOptimizer, func::Linear, set::T) where T <: LinSets
    cfunc = MOIU.canonical(func)
    addlinearconstraint!(m, cfunc, set)
    m.last_constraint_reference += 1
    ref = LCI{T}(m.last_constraint_reference)
    dict = constrdict(m, ref)
    dict[ref] = get_number_linear_constraints(m)
    push!(m.constraint_primal_solution, NaN)
    push!(m.constraint_dual_solution, NaN)
    push!(m.constraint_constant, func.constant)
    return ref
end

function addlinearconstraint!(m::LinQuadOptimizer, func::Linear, set::S) where S <: Union{LE, GE, EQ}
    addlinearconstraint!(m, func, lqs_char(m,set), _getrhs(set))
end

function addlinearconstraint!(m::LinQuadOptimizer, func::Linear, set::IV)
    addlinearconstraint!(m, func, lqs_char(m,set), set.lower)
    lqs_chgrngval!(m, [get_number_linear_constraints(m)], [set.upper - set.lower])
end

function addlinearconstraint!(m::LinQuadOptimizer, func::Linear, sense::Cchar, rhs)
    if abs(func.constant) > eps(Float64)
        warn("Constant in scalar function moved into set.")
    end
    add_linear_constraints!(m, [1], getcol.(m, func.variables), func.coefficients, [sense], [rhs - func.constant])
end

#=
    Add linear constraints (plural)
=#

function MOI.addconstraints!(m::LinQuadOptimizer, func::Vector{Linear}, set::Vector{S}) where S <: LinSets
    # canonicalize
    cfunc = MOIU.canonical.(func)

    @assert length(cfunc) == length(set)
    numrows = get_number_linear_constraints(m)
    addlinearconstraints!(m, cfunc, set)
    crefs = Vector{LCI{S}}(length(cfunc))
    for i in 1:length(cfunc)
        m.last_constraint_reference += 1
        ref = LCI{S}(m.last_constraint_reference)
        dict = constrdict(m, ref)
        dict[ref] = numrows + i
        push!(m.constraint_primal_solution, NaN)
        push!(m.constraint_dual_solution, NaN)
        push!(m.constraint_constant, cfunc[i].constant)
        crefs[i] = ref
    end
    return crefs
end

function addlinearconstraints!(m::LinQuadOptimizer, func::Vector{Linear}, set::Vector{S}) where S <: LinSets
    addlinearconstraints!(m, func, lqs_char.(m,set), [_getrhs(s) for s in set])
end

function addlinearconstraints!(m::LinQuadOptimizer, func::Vector{Linear}, set::Vector{IV})
    numrows = get_number_linear_constraints(m)
    addlinearconstraints!(m, func, lqs_char.(m,set), [s.lower for s in set])
    numrows2 = get_number_linear_constraints(m)
    lqs_chgrngval!(m, collect(numrows+1:numrows2), [s.upper - s.lower for s in set])
end

function addlinearconstraints!(m::LinQuadOptimizer, func::Vector{Linear}, sense::Vector{Cchar}, rhs::Vector{Float64})
    # loop through once to get number of non-zeros and to move rhs across
    nnz = 0
    for (i, f) in enumerate(func)
        if abs(f.constant) > eps(Float64)
            warn("Constant in scalar function moved into set.")
            rhs[i] -= f.constant
        end
        nnz += length(f.coefficients)
    end
    rowbegins = Vector{Int}(length(func))   # index of start of each row
    column_indices = Vector{Int}(nnz)       # flattened columns for each function
    nnz_vals = Vector{Float64}(nnz)         # corresponding non-zeros
    cnt = 1
    for (fi, f) in enumerate(func)
        rowbegins[fi] = cnt
        for (var, coef) in zip(f.variables, f.coefficients)
            column_indices[cnt] = getcol(m, var)
            nnz_vals[cnt] = coef
            cnt += 1
        end
    end
    add_linear_constraints!(m, rowbegins, column_indices, nnz_vals, sense, rhs)
end

#=
    Constraint set of Linear function
=#

MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintSet, ::Type{LCI{S}}) where S <: Union{LE, GE, EQ} = true
function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintSet, c::LCI{S}) where S <: Union{LE, GE, EQ}
    rhs = get_rhs(m, m[c])
    S(rhs+m.constraint_constant[m[c]])
end

#=
    Constraint function of Linear function
=#

MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintFunction, ::Type{<:LCI{<: LinSets}}) = true
function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintFunction, c::LCI{<: LinSets})
    # TODO more efficiently
    colidx, coefs = get_linear_constraint(m, m[c])
    Linear(m.variable_references[colidx+1], coefs, -m.constraint_constant[m[c]])
end

#=
    Scalar Coefficient Change of Linear Constraint
=#

MOI.canmodifyconstraint(m::LinQuadOptimizer, c::LCI{<: LinSets}, ::Type{MOI.ScalarCoefficientChange{Float64}}) = true
function MOI.modifyconstraint!(m::LinQuadOptimizer, c::LCI{<: LinSets}, chg::MOI.ScalarCoefficientChange{Float64})
    col = m.variable_mapping[chg.variable]
    lqs_chgcoef!(m, m[c], col, chg.new_coefficient)
end

#=
    Change RHS of linear constraint without modifying sense
=#

MOI.canmodifyconstraint(m::LinQuadOptimizer, c::LCI{S}, ::Type{S}) where S <: Union{LE, GE, EQ} = true
function MOI.modifyconstraint!(m::LinQuadOptimizer, c::LCI{S}, newset::S) where S <: Union{LE, GE, EQ}
    # the column 0 (or -1 in 0-index) is the rhs.
    lqs_chgcoef!(m, m[c], 0, _getrhs(newset))
end

MOI.canmodifyconstraint(m::LinQuadOptimizer, c::LCI{IV}, ::Type{IV}) = true
function MOI.modifyconstraint!(m::LinQuadOptimizer, c::LCI{IV}, set::IV)
    # the column 0 (or -1 in 0-index) is the rhs.
    # a range constraint has the RHS value of the lower limit of the range, and
    # a rngval equal to upper-lower.
    row = m[c]
    lqs_chgcoef!(m, row, 0, set.lower)
    lqs_chgrngval!(m, [row], [set.upper - set.lower])
end

#=
    Delete a linear constraint
=#

MOI.candelete(m::LinQuadOptimizer, c::LCI{<: LinSets}) = true
function MOI.delete!(m::LinQuadOptimizer, c::LCI{<: LinSets})
    deleteconstraintname!(m, c)
    dict = constrdict(m, c)
    row = dict[c]
    lqs_delrows!(m, row, row)
    deleteat!(m.constraint_primal_solution, row)
    deleteat!(m.constraint_dual_solution, row)
    deleteat!(m.constraint_constant, row)
    deleteref!(m, row, c)
end
function deleteref!(m::LinQuadOptimizer, row::Int, ref::LCI{<: LinSets})
    deleteref!(cmap(m).less_than, row, ref)
    deleteref!(cmap(m).greater_than, row, ref)
    deleteref!(cmap(m).equal_to, row, ref)
    deleteref!(cmap(m).interval, row, ref)
end

#=
    Transform scalar constraint
=#

function MOI.cantransformconstraint(m::LinQuadOptimizer, ref::LCI{S}, newset::S) where S
    false
end
function MOI.transformconstraint!(m::LinQuadOptimizer, ref::LCI{S}, newset::S) where S
    error("Cannot transform constraint of same set. use `modifyconstraint!` instead.")
end

function MOI.cantransformconstraint(m::LinQuadOptimizer, ref::LCI{S1}, newset::S2) where S1 where S2 <: Union{LE, GE, EQ}
    true
end
function MOI.transformconstraint!(m::LinQuadOptimizer, ref::LCI{S1}, newset::S2) where S1 where S2 <: Union{LE, GE, EQ}
    dict = constrdict(m, ref)
    row = dict[ref]
    lqs_chgsense!(m, [row], [lqs_char(m,newset)])
    m.last_constraint_reference += 1
    ref2 = LCI{S2}(m.last_constraint_reference)
    dict2 = constrdict(m, ref2)
    dict2[ref2] = row
    delete!(dict, ref)
    deleteconstraintname!(m, ref)
    return ref2
end
