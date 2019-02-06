#=
    Helper functions
=#

"""
    get_column(model::LinQuadOptimizer, index::MOI.VariableIndex)

Return the column of the variable `index` in `model`.
"""
function get_column(model::LinQuadOptimizer, index::MOI.VariableIndex)
    return variable_cache(model, index).column
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
    return variable_cache(model, index).name
end

# The rules for MOI are:
# - Allow setting duplicate variable names
# - Throw an error on get if the name is a duplicate
# So we need to store two things.
# 1. a mapping from VariableIndex -> Name
# 2. a reverse mapping from Name -> set of VariableIndex's with that name
function MOI.set(model::LinQuadOptimizer, ::MOI.VariableName,
                 index::MOI.VariableIndex, name::String)
    cache = variable_cache(model, index)
    current_name = cache.name
    cache.name = name
    if current_name != ""
        # Remove `index` from the set of current name.
        pop!(model.variable_names_rev[current_name], index)
    end
    if name != ""
        # Name is non-default, so store the reverse look-up.
        if !haskey(model.variable_names_rev, name)
            model.variable_names_rev[name] = Set{MOI.VariableIndex}()
        end
        push!(model.variable_names_rev[name], index)
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
    return
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
    model.last_variable_reference += 1
    index = MOI.VariableIndex(model.last_variable_reference)
    model.variable_mapping[index] = VariableCache(
        column = MOI.get(model, MOI.NumberOfVariables()),
        lower = -Inf,
        upper = Inf,
        variable_type = CONTINUOUS,
        bound_type = FREE,
        name = ""
    )
    push!(model.variable_references, index)
    push!(model.variable_primal_solution, NaN)
    push!(model.variable_dual_solution, NaN)
    return index
end

#=
    Add multiple variables
=#

function MOI.add_variables(model::LinQuadOptimizer, number_to_add::Int)
    previous_vars = MOI.get(model, MOI.NumberOfVariables())
    add_variables!(model, number_to_add)
    variable_indices = MOI.VariableIndex[]
    sizehint!(variable_indices, number_to_add)
    for i in 1:number_to_add
        model.last_variable_reference += 1
        index = MOI.VariableIndex(model.last_variable_reference)
        push!(variable_indices, index)
        model.variable_mapping[index] = VariableCache(
            column = previous_vars + i,
            lower = -Inf,
            upper = Inf,
            variable_type = CONTINUOUS,
            bound_type = FREE,
            name = ""
        )
        push!(model.variable_references, index)
        push!(model.variable_primal_solution, NaN)
        push!(model.variable_dual_solution, NaN)
    end
    return variable_indices
end

#=
    Check if reference is valid
=#

function MOI.is_valid(model::LinQuadOptimizer, index::MOI.VariableIndex)
    return haskey(model.variable_mapping, index)
end

#=
    Delete a variable
=#

function MOI.delete(model::LinQuadOptimizer, index::MOI.VariableIndex)
    __assert_valid__(model, index)
    cache = variable_cache(model, index)
    column = cache.column
    # Delete variables in internal model.
    delete_variables!(model, column, column)
    # Delete solution cache.
    deleteat!(model.variable_primal_solution, column)
    deleteat!(model.variable_dual_solution, column)
    # Delete from ordered list of variables.
    deleteat!(model.variable_references, column)
    # Delete from reverse name dictionary.
    delete!(model.variable_names_rev, cache.name)
    # Delete from variable_mapping.
    delete!(model.variable_mapping, index)
    # Decrement the column index of all other variables.
    for (variable, cache_2) in model.variable_mapping
        if cache_2.column > column
            cache_2.column -= 1
        end
    end
    return
end

#=
    MIP starts
=#

function MOI.supports(::LinQuadOptimizer, ::MOI.VariablePrimalStart,
                      ::Type{MOI.VariableIndex})
    return false
end

function MOI.set(model::LinQuadOptimizer, ::MOI.VariablePrimalStart,
                 indices::Vector{MOI.VariableIndex}, values::Vector{Float64})
    if !MOI.supports(model, MOI.VariablePrimalStart(), MOI.VariableIndex)
        throw(MOI.UnsupportedAttribute(MOI.VariablePrimalStart()))
    end
    add_mip_starts!(model, get_column.(Ref(model), indices), values)
    return
end

function MOI.set(model::LinQuadOptimizer, ::MOI.VariablePrimalStart,
                 index::MOI.VariableIndex, value::Float64)
    MOI.set(model, MOI.VariablePrimalStart(), [index], [value])
end
