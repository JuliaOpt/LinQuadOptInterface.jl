#=
    Vector valued constraints
=#
constrdict(model::LinQuadOptimizer, ::VLCI{MOI.Nonnegatives}) = cmap(model).nonnegatives
constrdict(model::LinQuadOptimizer, ::VLCI{MOI.Nonpositives}) = cmap(model).nonpositives
constrdict(model::LinQuadOptimizer, ::VLCI{MOI.Zeros})        = cmap(model).zeros

function MOI.addconstraint!(model::LinQuadOptimizer, func::VecLin, set::S) where S <: VecLinSets
    __assert_supported_constraint__(model, VecLin, S)
    @assert MOI.dimension(set) == length(func.constants)
    nrows = get_number_linear_constraints(model)
    add_linear_constraint(model, func, backend_type(model, set))
    nrows2 = get_number_linear_constraints(model)
    model.last_constraint_reference += 1
    index = VLCI{S}(model.last_constraint_reference)
    dict = constrdict(model, index)
    dict[index] = collect(nrows+1:nrows2)
    for i in 1:MOI.dimension(set)
        push!(model.constraint_primal_solution, NaN)
        push!(model.constraint_dual_solution, NaN)
        push!(model.constraint_constant, func.constants[i])
    end
    return index
end

function add_linear_constraint(model::LinQuadOptimizer, func::VecLin, sense::Cchar)

    # sort terms by output index, then by variable index
    func = copy(func)
    MOIU.sort_and_compress!(func.terms, MOIU.termindices, t -> true, MOIU.unsafe_add)
    columns       = [get_column(model, term.scalar_term.variable_index) for term in func.terms]
    coefficients  = [term.scalar_term.coefficient for term in func.terms]


    # Compute the row pointers in the compressed sparse row matrix. The row
    # pointers are defined recursively:
    #  r[1] = 1
    #  r[i] = r[i - 1] + (number of nonzero elements in the (i - 1)th row)
    #
    # To compute this, we first count up the number of nonzero elements in each
    # row (i - 1), storing the result in r[i]. Then we perform a cumsum on r,
    # storing the result back in r.
    num_rows = length(func.constants)
    row_pointers = fill(0, num_rows)
    row_pointers[1] = 1
    for term in func.terms
        row = term.output_index
        if row < 1 || row > num_rows
            throw(ArgumentError("Output index $row out of range [1, $num_rows]"))
        end
        if row < num_rows
            row_pointers[row + 1] += 1
        end
    end
    cumsum!(row_pointers, row_pointers)

    A = CSRMatrix{Float64}(row_pointers, columns, coefficients)
    add_linear_constraints!(model, A, fill(sense, num_rows), -func.constants)
end

function MOI.modify!(model::LinQuadOptimizer, index::VLCI{<: VecLinSets},
                     change::MOI.VectorConstantChange{Float64})
    rows = model[index]
    @assert length(change.new_constant) == length(rows)
    for (row, constant) in zip(rows, change.new_constant)
        change_rhs_coefficient!(model, row, -constant)
        model.constraint_constant[row] = constant
    end
end

function MOI.modify!(model::LinQuadOptimizer, index::VLCI{<: VecLinSets},
                     change::MOI.MultirowChange{Float64})
    column = get_column(model, change.variable)
    for (row, coef) in change.new_coefficients
        change_matrix_coefficient!(model, row, column, coef)
    end
end

function MOI.delete!(model::LinQuadOptimizer, index::VLCI{<:VecLinSets})
    __assert_valid__(model, index)
    delete_constraint_name(model, index)
    dict = constrdict(model, index)
    # we delete rows from largest to smallest here so that we don't have
    # to worry about updating references in a greater numbered row, only to
    # modify it later.
    for row in sort(dict[index], rev=true)
        delete_linear_constraints!(model, row, row)
        deleteat!(model.constraint_primal_solution, row)
        deleteat!(model.constraint_dual_solution, row)
        deleteat!(model.constraint_constant, row)
        # shift all the other references
        shift_references_after_delete_affine!(model, row)
    end
    delete!(dict, index)
end


#=
    Constraint set of Linear function
=#

MOI.canget(::LinQuadOptimizer, ::MOI.ConstraintSet, ::Type{VLCI{S}}) where S <: VecLinSets = true
function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintSet, index::VLCI{S}) where S <: VecLinSets
    constraint_indices = model[index]
    S(length(constraint_indices))
end

#=
    Constraint function of Linear function
=#

MOI.canget(::LinQuadOptimizer, ::MOI.ConstraintFunction, ::Type{VLCI{S}}) where S <: VecLinSets = true
function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintFunction, index::VLCI{<: VecLinSets})
    rows = model[index]
    constants = Float64[]
    terms = MOI.VectorAffineTerm{Float64}[]
    for (i, row) in enumerate(rows)
        rhs = get_rhs(model, row)
        push!(constants, -rhs)
        columns, coefficients = get_linear_constraint(model, row)
        for (column, coefficient) in zip(columns, coefficients)
            push!(terms, MOI.VectorAffineTerm{Float64}(
                    i,
                    MOI.ScalarAffineTerm{Float64}(
                        coefficient, model.variable_references[column]
                    )
                )
            )
        end
    end
    return MOI.VectorAffineFunction(terms, constants)
end
