#=
    Helper functions
=#
getcol(m::LinQuadOptimizer, ref::VarInd) = m.variable_mapping[ref]
getcol(m::LinQuadOptimizer, v::SinVar) = getcol(m, v.variable)

#=
    Get variable names
=#
MOI.supports(model::LinQuadOptimizer, ::MOI.VariableName, ::Type{VarInd}) = true
MOI.canget(m::LinQuadOptimizer, ::MOI.VariableName, ::Type{VarInd}) = true
function MOI.get(m::LinQuadOptimizer, ::MOI.VariableName, ref::VarInd)
    if haskey(m.variable_names, ref)
        m.variable_names[ref]
    else
        ""
    end
end

function MOI.set!(m::LinQuadOptimizer, ::MOI.VariableName, ref::VarInd, name::String)
    if haskey(m.variable_names_rev, name)
        if m.variable_names_rev[name] != ref
            error("Duplicate variable name: $(name)")
        end
    elseif name != ""
        if haskey(m.variable_names, ref)
            # we're renaming an existing variable
            old_name = m.variable_names[ref]
            delete!(m.variable_names_rev, old_name)
        end
        m.variable_names[ref] = name
        m.variable_names_rev[name] = ref
    end
end

#=
    Get variable by name
=#
function MOI.canget(m::LinQuadOptimizer, ::Type{MOI.VariableIndex}, name::String)
    haskey(m.variable_names_rev, name)
end
function MOI.get(m::LinQuadOptimizer, ::Type{MOI.VariableIndex}, name::String)
    return m.variable_names_rev[name]
end

#=
    Get number of variables
=#
MOI.canget(m::LinQuadOptimizer, ::MOI.NumberOfVariables) = true
function MOI.get(m::LinQuadOptimizer, ::MOI.NumberOfVariables)
    # TODO(odow): just return length(m.variable_references)?
    return get_number_variables(m)
end

#=
    List of Variable References
=#
MOI.canget(m::LinQuadOptimizer, ::MOI.ListOfVariableIndices) = true
function MOI.get(m::LinQuadOptimizer, ::MOI.ListOfVariableIndices)
    return m.variable_references
end

#=
    Add a single variable
=#

function MOI.addvariable!(m::LinQuadOptimizer)
    add_variables!(m, 1)
    # assumes we add columns linearly
    m.last_variable_reference += 1
    ref = VarInd(m.last_variable_reference)
    m.variable_mapping[ref] = MOI.get(m, MOI.NumberOfVariables())
    push!(m.variable_references, ref)
    push!(m.variable_primal_solution, NaN)
    push!(m.variable_dual_solution, NaN)
    return ref
end

#=
    Add multiple variables
=#

function MOI.addvariables!(m::LinQuadOptimizer, n::Int)
    previous_vars = MOI.get(m, MOI.NumberOfVariables())
    add_variables!(m, n)
    variable_references = VarInd[]
    sizehint!(variable_references, n)
    for i in 1:n
        # assumes we add columns linearly
        m.last_variable_reference += 1
        ref = VarInd(m.last_variable_reference)
        push!(variable_references, ref)
        # m.last_variable_reference might not be the column if we have
        # deleted variables.
        m.variable_mapping[ref] = previous_vars + i
        push!(m.variable_references, ref)
        push!(m.variable_primal_solution, NaN)
        push!(m.variable_dual_solution, NaN)
    end
    return variable_references
end

#=
    Check if reference is valid
=#

function MOI.isvalid(m::LinQuadOptimizer, ref::VarInd)
    if haskey(m.variable_mapping, ref)
        column = m.variable_mapping[ref]
        if column > 0 && column <= MOI.get(m, MOI.NumberOfVariables())
            return true
        end
    end
    return false
end

#=
    Delete a variable
=#

function MOI.delete!(m::LinQuadOptimizer, ref::VarInd)
    if !MOI.isvalid(m, ref)
        throw(MOI.InvalidIndex(ref))
    end
    col = m.variable_mapping[ref]
    delete_variables!(m, col, col)

    # delete from problem
    deleteat!(m.variable_references, col)
    deleteat!(m.variable_primal_solution, col)
    deleteat!(m.variable_dual_solution, col)
    deleteref!(m.variable_mapping, col, ref)

    # delete name if not ""
    if haskey(m.variable_names, ref)
        name = m.variable_names[ref]
        delete!(m.variable_names_rev, name)
        delete!(m.variable_names, ref)
    end

    # delete any bounds
    # deleting from a dict without the key does nothing
    deletebyval!(cmap(m).upper_bound, ref)
    deletebyval!(cmap(m).lower_bound, ref)
    deletebyval!(cmap(m).fixed_bound, ref)
    deletebyval!(cmap(m).interval_bound, ref)

end

# temp fix - change storage for bounds
# TODO(@joaquimg): what did you mean by this?
function deletebyval!(dict::Dict{S,T}, index::T) where {S,T}
    for (key, val) in dict
        if val == index
            delete!(dict, key)
        end
    end
end

#=
    MIP starts
=#
MOI.supports(model::LinQuadOptimizer, ::MOI.VariablePrimalStart, ::Type{VarInd}) = false

function MOI.set!(m::LinQuadOptimizer, ::MOI.VariablePrimalStart, ref::VarInd, val::Float64)
    add_mip_starts!(m, [getcol(m, ref)], [val])
end

function MOI.set!(m::LinQuadOptimizer, ::MOI.VariablePrimalStart, refs::Vector{VarInd}, vals::Vector{Float64})
    add_mip_starts!(m, getcol.(Ref(m), refs), vals)
end
