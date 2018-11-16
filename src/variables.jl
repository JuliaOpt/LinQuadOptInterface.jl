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
function MOI.get(model::LinQuadOptimizer, ::MOI.VariableName, index::VarInd)
    if haskey(model.variable_names, index)
        return model.variable_names[index]
    else
        return ""
    end
end

# The rules for MOI are:
# - Allow setting duplicate variable names
# - Throw an error on get if the name is a duplicate
# So we need to store two things.
# 1. a mapping from VariableIndex -> Name
# 2. a reverse mapping from Name -> set of VariableIndex's with that name
function MOI.set(model::LinQuadOptimizer,
                 ::MOI.VariableName, index::VarInd, name::String)
    if haskey(model.variable_names, index)
        # This variable already has a name, we must be changing it.
        current_name = model.variable_names[index]
        # Remove `index` from the set of current name.
        pop!(model.variable_names_rev[current_name], index)
    end
    if name != ""
        # We're changing the name to something non-default, so store it.
        model.variable_names[index] = name
        if !haskey(model.variable_names_rev, name)
            model.variable_names_rev[name] = Set{VarInd}()
        end
        push!(model.variable_names_rev[name], index)
    else
        # We're changing the name to the default, so we don't store it. Note
        # that if `model.variable_names` doesn't have a key `index`, this does
        # nothing.
        delete!(model.variable_names, index)
    end
    return
end

#=
    Get variable by name
=#

function MOI.get(model::LinQuadOptimizer, ::Type{MOI.VariableIndex}, name::String)
    if haskey(model.variable_names_rev, name)
        variable_set = model.variable_names_rev[name]
        if length(variable_set) == 1
            return first(model.variable_names_rev[name])
        elseif length(variable_set) > 1
            error("Cannot get variable because the name $(name) is a duplicate.")
        end
    end
    return nothing
end

#=
    Get number of variables
=#

function MOI.get(model::LinQuadOptimizer, ::MOI.NumberOfVariables)
    # TODO(odow): work out why this is not the same as
    # length(model.variable_references)
    return get_number_variables(model)
end

#=
    List of Variable References
=#

function MOI.get(model::LinQuadOptimizer, ::MOI.ListOfVariableIndices)
    return model.variable_references
end

#=
    Add a single variable
=#

function MOI.add_variable(model::LinQuadOptimizer)
    add_variables!(model, 1)
    # assumes we add columns linearly
    model.last_variable_reference += 1
    index = VarInd(model.last_variable_reference)
    model.variable_mapping[index] = MOI.get(model, MOI.NumberOfVariables())
    push!(model.variable_references, index)
    push!(model.variable_primal_solution, NaN)
    push!(model.variable_dual_solution, NaN)
    model.variable_type[index] = Continuous
    return index
end

#=
    Add multiple variables
=#

function MOI.add_variables(model::LinQuadOptimizer, number_to_add::Int)
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
        model.variable_type[index] = Continuous
    end
    return variable_indices
end

#=
    Check if reference is valid
=#

function MOI.is_valid(model::LinQuadOptimizer, index::VarInd)
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

function MOI.delete(model::LinQuadOptimizer, index::VarInd)
    __assert_valid__(model, index)
    column = get_column(model, index)
    delete_variables!(model, column, column)
    # delete from problem
    deleteat!(model.variable_references, column)
    delete!(model.variable_type, index)
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
function MOI.supports(
        ::LinQuadOptimizer, ::MOI.VariablePrimalStart, ::Type{VarInd})
    return false
end

function MOI.set(model::LinQuadOptimizer, ::MOI.VariablePrimalStart,
                  indices::Vector{VarInd}, values::Vector{Float64})
    if !MOI.supports(model, MOI.VariablePrimalStart(), MOI.VariableIndex)
        throw(MOI.UnsupportedAttribute(MOI.VariablePrimalStart()))
    end
    add_mip_starts!(model, get_column.(Ref(model), indices), values)
end

function MOI.set(model::LinQuadOptimizer, ::MOI.VariablePrimalStart,
                  index::VarInd, value::Float64)
    MOI.set(model, MOI.VariablePrimalStart(), [index], [value])
end
