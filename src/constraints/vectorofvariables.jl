#=
    Vector valued bounds
        VectorOfVariables -in- Nonnegatives
        VectorOfVariables -in- Nonpositives
        VectorOfVariables -in- Zeros

    SOS
        VectorOfVariables -in- SOS1
        VectorOfVariables -in- SOS2
=#
constrdict(model::LinQuadOptimizer, ::VVCI{MOI.Nonnegatives}) = cmap(model).vv_nonnegatives
constrdict(model::LinQuadOptimizer, ::VVCI{MOI.Nonpositives}) = cmap(model).vv_nonpositives
constrdict(model::LinQuadOptimizer, ::VVCI{MOI.Zeros}) = cmap(model).vv_zeros

constrdict(model::LinQuadOptimizer, ::VVCI{SOS1}) = cmap(model).sos1
constrdict(model::LinQuadOptimizer, ::VVCI{SOS2}) = cmap(model).sos2

function MOI.add_constraint(model::LinQuadOptimizer, func::VecVar, set::S) where S <: VecLinSets
    __assert_supported_constraint__(model, VecVar, S)
    @assert length(func.variables) == MOI.dimension(set)
    rows = get_number_linear_constraints(model)
    num_vars = MOI.dimension(set)
    add_linear_constraints!(model,
        CSRMatrix{Float64}(collect(1:num_vars+1),
                           get_column.(Ref(model), func.variables),
                           ones(num_vars)),
        fill(backend_type(model, set), num_vars),
        zeros(num_vars)
    )
    model.last_constraint_reference += 1
    index = VVCI{S}(model.last_constraint_reference)
    dict = constrdict(model, index)
    dict[index] = collect((rows + 1):(rows + num_vars))
    append!(model.constraint_primal_solution, fill(NaN, num_vars))
    append!(model.constraint_dual_solution, fill(NaN, num_vars))
    append!(model.constraint_constant, fill(0.0, num_vars))
    return index
end

function MOI.delete(model::LinQuadOptimizer, index::VVCI{S}) where S <: VecLinSets
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
    Get constraint set of vector variable bound
=#

function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintSet, index::VVCI{S}) where S <: VecLinSets
    S(length(model[index]))
end

#=
    Get constraint function of vector variable bound (linear ctr)
=#

function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintFunction, index::VVCI{<: VecLinSets})
    rows = model[index]
    variables = VarInd[]
    sizehint!(variables, length(rows))
    for row in rows
        cols, coefs = get_linear_constraint(model, row)
        if length(cols) != 1
            error("Unexpected constraint")
        end
        push!(variables, model.variable_references[cols[1]])
    end
    return VecVar(variables)
end

#=
    SOS constraints
=#

function MOI.add_constraint(model::LinQuadOptimizer, func::VecVar, set::S) where S <: Union{MOI.SOS1, MOI.SOS2}
    __assert_supported_constraint__(model, VecVar, S)
    make_problem_type_integer(model)
    add_sos_constraint!(model, get_column.(Ref(model), func.variables),
                        set.weights, backend_type(model, set))
    model.last_constraint_reference += 1
    index = VVCI{S}(model.last_constraint_reference)
    dict = constrdict(model, index)
    dict[index] = length(cmap(model).sos1) + length(cmap(model).sos2) + 1
    return index
end

function MOI.delete(model::LinQuadOptimizer, index::VVCI{<:Union{SOS1, SOS2}})
    __assert_valid__(model, index)
    delete_constraint_name(model, index)
    dict = constrdict(model, index)
    sos_row = dict[index]
    delete_sos!(model, sos_row, sos_row)
    deleteref!(cmap(model).sos1, sos_row, index)
    deleteref!(cmap(model).sos2, sos_row, index)
    if !has_integer(model)
        make_problem_type_continuous(model)
    end
end

function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintSet, index::VVCI{S}) where S <: Union{MOI.SOS1, MOI.SOS2}
    indices, weights, types = get_sos_constraint(model, model[index])
    set = S(weights)
    @assert types == backend_type(model, set)
    return set
end

function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintFunction, index::VVCI{<:Union{SOS1, SOS2}})
    indices, weights, types = get_sos_constraint(model, model[index])
    return VecVar(model.variable_references[indices])
end
