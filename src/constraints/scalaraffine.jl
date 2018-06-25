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
    addlinearconstraint!(m, func, backend_type(m,set), _getrhs(set))
end

function addlinearconstraint!(m::LinQuadOptimizer, func::Linear, set::IV)
    columns = [getcol(m, term.variable_index) for term in func.terms]
    coefficients = [term.coefficient for term in func.terms]
    A = CSRMatrix{Float64}([1], columns, coefficients)
    add_ranged_constraints!(m, A, [set.lower], [set.upper])
end

function addlinearconstraint!(m::LinQuadOptimizer, func::Linear, sense::Cchar, rhs)
    if abs(func.constant) > eps(Float64)
        warn("Constant in scalar function moved into set.")
    end
    columns = [getcol(m, term.variable_index) for term in func.terms]
    coefficients = [term.coefficient for term in func.terms]
    A = CSRMatrix{Float64}([1], columns, coefficients)
    add_linear_constraints!(m, A, [sense], [rhs - func.constant])
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
    addlinearconstraints!(m, func, backend_type.(m,set), [_getrhs(s) for s in set])
end

function addlinearconstraints!(m::LinQuadOptimizer, func::Vector{Linear}, set::Vector{IV})
    # loop through once to get number of non-zeros and to move rhs across
    lowerbounds = [s.lower for s in set]
    upperbounds = [s.upper for s in set]
    nnz = 0
    for (i, f) in enumerate(func)
        if abs(f.constant) > eps(Float64)
            warn("Constant in scalar function moved into set.")
            lowerbounds[i] -= f.constant
            upperbounds[i] -= f.constant
        end
        nnz += length(f.terms)
    end
    row_pointers = Vector{Int}(length(func))  # index of start of each row
    columns = Vector{Int}(nnz)                # flattened columns for each function
    coefficients  = Vector{Float64}(nnz)      # corresponding non-zeros
    i = 1
    for (fi, f) in enumerate(func)
        row_pointers[fi] = i
        for term in f.terms
            columns[i] = getcol(m, term.variable_index)
            coefficients[i] = term.coefficient
            i += 1
        end
    end
    A = CSRMatrix{Float64}(row_pointers, columns, coefficients)
    add_ranged_constraints!(m, A, lowerbounds, upperbounds)
end

function addlinearconstraints!(m::LinQuadOptimizer, func::Vector{Linear}, sense::Vector{Cchar}, rhs::Vector{Float64})
    # loop through once to get number of non-zeros and to move rhs across
    nnz = 0
    for (i, f) in enumerate(func)
        if abs(f.constant) > eps(Float64)
            warn("Constant in scalar function moved into set.")
            rhs[i] -= f.constant
        end
        nnz += length(f.terms)
    end
    row_pointers = Vector{Int}(length(func))  # index of start of each row
    columns = Vector{Int}(nnz)                # flattened columns for each function
    coefficients = Vector{Float64}(nnz)       # corresponding non-zeros
    i = 1
    for (row, f) in enumerate(func)
        row_pointers[row] = i
        for term in f.terms
            columns[i] = getcol(m, term.variable_index)
            coefficients[i] = term.coefficient
            i += 1
        end
    end
    A = CSRMatrix{Float64}(row_pointers, columns, coefficients)
    add_linear_constraints!(m, A, sense, rhs)
end

#=
    Constraint set of Linear function
=#

MOI.canget(::LinQuadOptimizer, ::MOI.ConstraintSet, ::Type{LCI{S}}) where S <: Union{LE, GE, EQ} = true
function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintSet, c::LCI{S}) where S <: Union{LE, GE, EQ}
    rhs = get_rhs(m, m[c])
    S(rhs+m.constraint_constant[m[c]])
end

MOI.canget(::LinQuadOptimizer, ::MOI.ConstraintSet, ::Type{LCI{IV}}) = true
function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintSet, c::LCI{IV})
    lowerbound, upperbound = get_range(m, m[c])
    IV(lowerbound+m.constraint_constant[m[c]], upperbound + m.constraint_constant[m[c]])
end

#=
    Constraint function of Linear function
=#

MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintFunction, ::Type{<:LCI{<: LinSets}}) = true
function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintFunction, c::LCI{<: LinSets})
    # TODO more efficiently
    colidx, coefs = get_linear_constraint(m, m[c])
    terms = map(
        (v,c)->MOI.ScalarAffineTerm{Float64}(c,v),
        m.variable_references[colidx],
        coefs
    )
    Linear(terms, -m.constraint_constant[m[c]])
end

#=
    Scalar Coefficient Change of Linear Constraint
=#

MOI.canmodify(m::LinQuadOptimizer, ::Type{LCI{S}}, ::Type{MOI.ScalarCoefficientChange{Float64}}) where S <: LinSets = true
function MOI.modify!(m::LinQuadOptimizer, c::LCI{S}, chg::MOI.ScalarCoefficientChange{Float64}) where S <: LinSets
    col = m.variable_mapping[chg.variable]
    change_matrix_coefficient!(m, m[c], col, chg.new_coefficient)
end

#=
    Change RHS of linear constraint without modifying sense
=#

MOI.canset(m::LinQuadOptimizer, ::MOI.ConstraintSet, ::Type{LCI{S}}) where S <: Union{LE, GE, EQ} = true
function MOI.set!(m::LinQuadOptimizer, ::MOI.ConstraintSet, c::LCI{S}, newset::S) where S <: Union{LE, GE, EQ}
    change_rhs_coefficient!(m, m[c], _getrhs(newset))
end

MOI.canset(m::LinQuadOptimizer, ::MOI.ConstraintSet, ::Type{LCI{IV}}) = true
function MOI.set!(m::LinQuadOptimizer, ::MOI.ConstraintSet, c::LCI{IV}, set::IV)
    modify_ranged_constraints!(m, [m[c]], [set.lower], [set.upper])
end

#=
    Delete a linear constraint
=#

MOI.candelete(m::LinQuadOptimizer, c::LCI{<: LinSets}) = MOI.isvalid(m, c)
function MOI.delete!(m::LinQuadOptimizer, c::LCI{<: LinSets})
    deleteconstraintname!(m, c)
    dict = constrdict(m, c)
    row = dict[c]
    delete_linear_constraints!(m, row, row)
    deleteat!(m.constraint_primal_solution, row)
    deleteat!(m.constraint_dual_solution, row)
    deleteat!(m.constraint_constant, row)
    # shift all the other references
    shift_references_after_delete_affine!(m, row)
    delete!(dict, c)
end

#=
    Transform scalar constraint
=#

function MOI.cantransform(m::LinQuadOptimizer, ref::LCI{S}, newset::S) where S
    false
end
function MOI.transform!(m::LinQuadOptimizer, ::LCI{S}, newset::S) where S
    error("Cannot transform constraint of same set. use `set!` instead.")
end

function MOI.cantransform(m::LinQuadOptimizer, ref::LCI{S1}, newset::S2) where S1 where S2 <: Union{LE, GE, EQ}
    true
end
function MOI.transform!(m::LinQuadOptimizer, ref::LCI{S1}, newset::S2) where S1 where S2 <: Union{LE, GE, EQ}
    dict = constrdict(m, ref)
    row = dict[ref]
    change_linear_constraint_sense!(m, [row], [backend_type(m,newset)])
    m.last_constraint_reference += 1
    ref2 = LCI{S2}(m.last_constraint_reference)
    dict2 = constrdict(m, ref2)
    dict2[ref2] = row
    delete!(dict, ref)
    deleteconstraintname!(m, ref)
    return ref2
end
