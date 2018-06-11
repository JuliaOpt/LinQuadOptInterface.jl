#=
        ScalarQuadraticFunction -in- LessThan
        ScalarQuadraticFunction -in- GreaterThan
        ScalarQuadraticFunction -in- EqualTo

    TODO:
        - constraint sets
        - constraint functions
        - deleting constraints
=#
constrdict(m::LinQuadOptimizer, ::QCI{LE})  = cmap(m).q_less_than
constrdict(m::LinQuadOptimizer, ::QCI{GE})  = cmap(m).q_greater_than
constrdict(m::LinQuadOptimizer, ::QCI{EQ})  = cmap(m).q_equal_to

function MOI.addconstraint!(m::LinQuadOptimizer, func::Quad, set::S) where S <: Union{LE, GE, EQ}
    addquadraticconstraint!(m, func, set)
    m.last_constraint_reference += 1
    ref = QCI{S}(m.last_constraint_reference)
    dict = constrdict(m, ref)
    push!(m.qconstraint_primal_solution, NaN)
    push!(m.qconstraint_dual_solution, NaN)
    # dict[ref] = get_number_quadratic_constraints(m)
    dict[ref] = length(m.qconstraint_primal_solution)
    return ref
end

function addquadraticconstraint!(m::LinQuadOptimizer, func::Quad, set::S) where S<: Union{LE, GE, EQ}
    addquadraticconstraint!(m, func, backend_type(m,set), _getrhs(set))
end

function addquadraticconstraint!(m::LinQuadOptimizer, f::Quad, sense::Cchar, rhs::Float64)
    if abs(f.constant) > 0
        warn("Constant in quadratic function. Moving into set")
    end
    quadratic_columns_1 = [getcol(m, term.variable_index_1) for term in f.quadratic_terms]
    quadratic_columns_2 = [getcol(m, term.variable_index_2) for term in f.quadratic_terms]
    quadratic_coefficients = [term.coefficient for term in f.quadratic_terms]
    ri, ci, vi = reduce_duplicates!(
        quadratic_columns_1,
        quadratic_columns_2,
        quadratic_coefficients
    )
    affine_columns = [getcol(m, term.variable_index) for term in f.affine_terms]
    affine_coefficients = [term.coefficient for term in f.affine_terms]
    add_quadratic_constraint!(m,
        affine_columns,
        affine_coefficients,
        rhs - f.constant,
        sense,
        ri, ci, vi
    )
end

"""
    reduce_duplicates!(rows::Vector{T}, cols::Vector{T}, vals::Vector{S})

Given a matrix specified by row indices in `rows`, column indices in `cols` and
coefficients in `vals`, return new `rows`, `cols`, and `vals` vectors with
duplicate elements summed and any coefficients in the lower triangle moved to
the upper triangle.

This function swaps element `i` in `rows` and `cols` if `rows[i]>cols[i]`.

# Examples
```jldoctest
julia> reduce_duplicates!(
    [1,   2, 2, 2],  # rows
    [1,   1, 2, 2],  # cols
    [1, 0.5, 1, 1]   # vals
    )
([1, 1, 2], [1, 2, 2], [1.0, 0.5, 2.0])
```
"""
function reduce_duplicates!(rows::Vector{T}, cols::Vector{T}, vals::Vector{S}) where T where S
    @assert length(rows) == length(cols) == length(vals)
    for i in 1:length(rows)
        if rows[i] > cols[i]
            tmp = rows[i]
            rows[i] = cols[i]
            cols[i] = tmp
        end
    end
    findnz(sparse(rows, cols, vals))
end


MOI.candelete(m::LinQuadOptimizer, c::QCI{<: LinSets}) = MOI.isvalid(m, c)
function MOI.delete!(m::LinQuadOptimizer, c::QCI{<: LinSets})
    deleteconstraintname!(m, c)
    dict = constrdict(m, c)
    row = dict[c]
    delete_quadratic_constraints!(m, row, row)
    deleteat!(m.qconstraint_primal_solution, row)
    deleteat!(m.qconstraint_dual_solution, row)
    # shift all the other references
    shift_references_after_delete_quadratic!(m, row)
    delete!(dict, c)
end

#=
    Constraint set of Linear function
=#

MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintSet, ::Type{QCI{S}}) where S <: Union{LE, GE, EQ} = true
function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintSet, c::QCI{S}) where S <: Union{LE, GE, EQ}
    rhs = get_rhs(m, m[c])
    S(rhs)
end

MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintSet, ::Type{QCI{IV}}) = false

#=
    Constraint function of Linear function
=#

MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintFunction, ::Type{<:QCI{<: LinSets}}) = true
function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintFunction, c::QCI{<: LinSets})
    # TODO more efficiently
    colidx, coefs, Q = get_quadratic_constraint(m, m[c])
    affine_terms = map(
        (v,c)->MOI.ScalarAffineTerm{Float64}(c,v),
        m.variable_references[colidx+1],
        coefs
    )
    rows = rowvals(Q)
    vals = nonzeros(Q)
    nrows, ncols = size(Q)
    quadratic_terms = MOI.ScalarQuadraticTerm{Float64}[]
    sizehint!(quadratic_terms, length(vals))
    for i = 1:ncols
        for j in nzrange(Q, i)
            row = rows[j]
            val = vals[j]
            push!(quadratic_terms, MOI.ScalarQuadraticTerm{Float64}(2*val, m.variable_references[row], m.variable_references[i]))
        end
    end
    Quad(affine_terms, quadratic_terms, 0.0) # constant was moved into set
end