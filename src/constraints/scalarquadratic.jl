#=
        ScalarQuadraticFunction -in- LessThan
        ScalarQuadraticFunction -in- GreaterThan
        ScalarQuadraticFunction -in- EqualTo

    TODO:
        - constraint sets
        - constraint functions
        - deleting constraints
=#
constrdict(model::LinQuadOptimizer, ::QCI{LE})  = cmap(model).q_less_than
constrdict(model::LinQuadOptimizer, ::QCI{GE})  = cmap(model).q_greater_than
constrdict(model::LinQuadOptimizer, ::QCI{EQ})  = cmap(model).q_equal_to

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

"""
    canonical_reduction(model::LinQuadOptimizer, func::Quad)

Reduce a ScalarQuadraticFunction into five arrays, returned in the following
order:
 1. a vector of affine column indices
 2. a vector of affine coefficients
 3. a vector of quadratic row indices
 4. a vector of quadratic column indices
 5. a vector of quadratic coefficients
"""
function canonical_reduction(model::LinQuadOptimizer, func::Quad)
    quad_columns_1, quad_columns_2, quad_coefficients = reduce_duplicates!(
        [get_column(model, term.variable_index_1) for term in func.quadratic_terms],
        [get_column(model, term.variable_index_2) for term in func.quadratic_terms],
        [term.coefficient for term in func.quadratic_terms]
    )
    affine_columns = [get_column(model, term.variable_index) for term in func.affine_terms]
    affine_coefficients = [term.coefficient for term in func.affine_terms]
    return affine_columns, affine_coefficients, quad_columns_1, quad_columns_2, quad_coefficients
end

function MOI.add_constraint(model::LinQuadOptimizer, func::Quad, set::S) where S <: Union{LE, GE, EQ}
    add_quadratic_constraint(model, func, set)
    model.last_constraint_reference += 1
    index = QCI{S}(model.last_constraint_reference)
    dict = constrdict(model, index)
    push!(model.qconstraint_primal_solution, NaN)
    push!(model.qconstraint_dual_solution, NaN)
    dict[index] = get_number_quadratic_constraints(model)
    return index
end

function add_quadratic_constraint(model::LinQuadOptimizer, func::Quad, set::S) where S<: Union{LE, GE, EQ}
    add_quadratic_constraint(model, func, backend_type(model, set),
                             MOIU.getconstant(set))
end

function add_quadratic_constraint(model::LinQuadOptimizer, func::Quad, sense, rhs::Float64)
    if abs(func.constant) > 0
        Compat.@warn("Constant in quadratic function. Moving into set")
    end
    (aff_cols, aff_coeffs, I, J, V) = canonical_reduction(model, func)
    add_quadratic_constraint!(model, aff_cols, aff_coeffs, rhs - func.constant,
                              sense, I, J, V)
end

function MOI.delete(model::LinQuadOptimizer, index::QCI{<: LinSets})
    __assert_valid__(model, index)
    delete_constraint_name(model, index)
    dict = constrdict(model, index)
    row = dict[index]
    delete_quadratic_constraints!(model, row, row)
    deleteat!(model.qconstraint_primal_solution, row)
    deleteat!(model.qconstraint_dual_solution, row)
    # shift all the other references
    shift_references_after_delete_quadratic!(model, row)
    delete!(dict, index)
end

#=
    Constraint set of Linear function
=#

function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintSet, index::QCI{S}) where S <: Union{LE, GE, EQ}
    rhs = get_quadratic_rhs(model, model[index])
    S(rhs)
end

#=
    Constraint function of Linear function
=#

function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintFunction, index::QCI{<: LinSets})
    aff_cols, aff_coeffs, Q = get_quadratic_constraint(model, model[index])
    affine_terms = map(
        (col, coef)->MOI.ScalarAffineTerm{Float64}(coef, model.variable_references[col]),
        aff_cols, aff_coeffs)
    rows = rowvals(Q)
    vals = nonzeros(Q)
    nrows, ncols = size(Q)
    quadratic_terms = MOI.ScalarQuadraticTerm{Float64}[]
    sizehint!(quadratic_terms, length(vals))
    for i = 1:ncols
        for j in nzrange(Q, i)
            row = rows[j]
            val = vals[j]
            push!(quadratic_terms, MOI.ScalarQuadraticTerm{Float64}(2*val,
                model.variable_references[row], model.variable_references[i]))
        end
    end
    Quad(affine_terms, quadratic_terms, 0.0) # constant was moved into set
end
