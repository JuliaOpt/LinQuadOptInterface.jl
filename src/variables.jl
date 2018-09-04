#=
    Helper functions
=#
"""
    get_column(model::LinQuadOptimizer, index::MOI.VariableIndex)

Return the column of the variable `index` in `model`.
"""
function get_column(model::LinQuadOptimizer, index::VarInd)
    return model.variable_mapping[index]
end

"""
    get_column(model::LinQuadOptimizer, variable::MOI.SingleVariable)

Return the column of the variable inside a `MOI.SingleVariable` function in
`model`.
"""
function get_column(model::LinQuadOptimizer, variable::SinVar)
    return get_column(model, variable.variable)
end

#=
    Get variable names
=#
MOI.supports(::LinQuadOptimizer, ::MOI.VariableName, ::Type{VarInd}) = true
MOI.canget(::LinQuadOptimizer, ::MOI.VariableName, ::Type{VarInd}) = true
function MOI.get(model::LinQuadOptimizer, ::MOI.VariableName, index::VarInd)
    if haskey(model.variable_names, index)
        return model.variable_names[index]
    else
        return ""
    end
end

function MOI.set!(model::LinQuadOptimizer, ::MOI.VariableName, index::VarInd, name::String)
    if haskey(model.variable_names_rev, name)
        if model.variable_names_rev[name] != index
            error("Duplicate variable name: $(name)")
        end
    elseif name != ""
        if haskey(model.variable_names, index)
            # we're renaming an existing variable
            old_name = model.variable_names[index]
            delete!(model.variable_names_rev, old_name)
        end
        model.variable_names[index] = name
        model.variable_names_rev[name] = index
    end
    return
end

#=
    Get variable by name
=#
function MOI.canget(model::LinQuadOptimizer, ::Type{MOI.VariableIndex}, name::String)
    return haskey(model.variable_names_rev, name)
end
function MOI.get(model::LinQuadOptimizer, ::Type{MOI.VariableIndex}, name::String)
    return model.variable_names_rev[name]
end

#=
    Get number of variables
=#
MOI.canget(::LinQuadOptimizer, ::MOI.NumberOfVariables) = true
function MOI.get(model::LinQuadOptimizer, ::MOI.NumberOfVariables)
    # TODO(odow): work out why this is not the same as
    # length(model.variable_references)
    return get_number_variables(model)
end

#=
    List of Variable References
=#
MOI.canget(::LinQuadOptimizer, ::MOI.ListOfVariableIndices) = true
function MOI.get(model::LinQuadOptimizer, ::MOI.ListOfVariableIndices)
    return model.variable_references
end

#=
    Add a single variable
=#

function MOI.addvariable!(model::LinQuadOptimizer)
    add_variables!(model, 1)
    # assumes we add columns linearly
    model.last_variable_reference += 1
    index = VarInd(model.last_variable_reference)
    model.variable_mapping[index] = MOI.get(model, MOI.NumberOfVariables())
    push!(model.variable_references, index)
    push!(model.variable_primal_solution, NaN)
    push!(model.variable_dual_solution, NaN)
    return index
end

#=
    Add multiple variables
=#

function MOI.addvariables!(model::LinQuadOptimizer, number_to_add::Int)
    previous_vars = MOI.get(model, MOI.NumberOfVariables())
    add_variables!(model, number_to_add)
    variable_indices = VarInd[]
    sizehint!(variable_indices, number_to_add)
    for i in 1:number_to_add
        # assumes we add columns linearly
        model.last_variable_reference += 1
        index = VarInd(model.last_variable_reference)
        push!(variable_indices, index)
        # m.last_variable_reference might not be the column if we have
        # deleted variables.
        model.variable_mapping[index] = previous_vars + i
        push!(model.variable_references, index)
        push!(model.variable_primal_solution, NaN)
        push!(model.variable_dual_solution, NaN)
    end
    return variable_indices
end

#=
    Check if reference is valid
=#

function MOI.isvalid(model::LinQuadOptimizer, index::VarInd)
    if haskey(model.variable_mapping, index)
        column = get_column(model, index)
        if 1 <= column <= MOI.get(model, MOI.NumberOfVariables())
            return true
        end
    end
    return false
end

#=
    Delete a variable
=#

function MOI.delete!(model::LinQuadOptimizer, index::VarInd)
    __assert_valid__(model, index)
    column = get_column(model, index)
    delete_variables!(model, column, column)
    # delete from problem
    deleteat!(model.variable_references, column)
    deleteat!(model.variable_primal_solution, column)
    deleteat!(model.variable_dual_solution, column)
    deleteref!(model.variable_mapping, column, index)
    # delete name if not ""
    if haskey(model.variable_names, index)
        name = model.variable_names[index]
        delete!(model.variable_names_rev, name)
        delete!(model.variable_names, index)
    end
    # Delete any bounds (deleting from a dict without the key does nothing).
    constraint_map = cmap(model)
    delete_from_dictionary_by_value(constraint_map.upper_bound, index)
    delete_from_dictionary_by_value(constraint_map.lower_bound, index)
    delete_from_dictionary_by_value(constraint_map.fixed_bound, index)
    delete_from_dictionary_by_value(constraint_map.interval_bound, index)
    return
end

"""
    delete_from_dictionary_by_value(dict::Dict{K, V}, index::V)

Delete the entries of a dictionary where the value of a `key=>value` pair is
equal to `index`.

This is used to remove bounds when deleting a variable.
"""
function delete_from_dictionary_by_value(dict::Dict{K, V}, index::V) where {K, V}
    for (key, value) in dict
        if value == index
            delete!(dict, key)
        end
    end
end

#=
    MIP starts
=#
function MOI.supports(::LinQuadOptimizer, ::MOI.VariablePrimalStart, ::Type{VarInd})
    return false
end

function MOI.set!(model::LinQuadOptimizer, ::MOI.VariablePrimalStart,
                  indices::Vector{VarInd}, values::Vector{Float64})
    add_mip_starts!(model, get_column.(Ref(model), indices), values)
end

function MOI.set!(model::LinQuadOptimizer, ::MOI.VariablePrimalStart,
                  index::VarInd, value::Float64)
    MOI.set!(model, MOI.VariablePrimalStart(), [index], [value])
end
