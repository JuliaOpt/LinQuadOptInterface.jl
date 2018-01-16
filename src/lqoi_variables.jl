#=
    Helper functions
=#
getcol(m::LinQuadSolverInstance, ref::VarInd) = m.variable_mapping[ref]
getcol(m::LinQuadSolverInstance, v::SinVar) = getcol(m, v.variable)

#=
    Get number of variables
=#

function MOI.get(m::LinQuadSolverInstance, ::MOI.NumberOfVariables)
    return lqs_getnumcols(m)
end
MOI.canget(m::LinQuadSolverInstance, ::MOI.NumberOfVariables) = true

#=
    Get variable references
=#

function MOI.get(m::LinQuadSolverInstance, ::MOI.ListOfVariableIndices)
    return m.variable_references
end
MOI.canget(m::LinQuadSolverInstance, ::MOI.ListOfVariableIndices) = true

#=
    Add a variable to the LinQuadSolverInstance
Returns a MathOptInterface VariableIndex. So we need to increment the
variable reference counter in the m.variablemapping, and store the column number
in the dictionary
=#

function MOI.addvariable!(m::LinQuadSolverInstance)
    lqs_newcols!(m, 1)
    # assumes we add columns linearly
    m.last_variable_reference += 1
    ref = MOI.VariableIndex(m.last_variable_reference)
    m.variable_mapping[ref] = MOI.get(m, MOI.NumberOfVariables())
    push!(m.variable_references, ref)
    push!(m.variable_primal_solution, NaN)
    push!(m.variable_dual_solution, NaN)
    return ref
end

function MOI.addvariables!(m::LinQuadSolverInstance, n::Int)
    previous_vars = MOI.get(m, MOI.NumberOfVariables())
    lqs_newcols!(m, n)
    variable_references = MOI.VariableIndex[]
    sizehint!(variable_references, n)
    for i in 1:n
        # assumes we add columns linearly
        m.last_variable_reference += 1
        ref = MOI.VariableIndex(m.last_variable_reference)
        push!(variable_references, ref)
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

function MOI.isvalid(m::LinQuadSolverInstance, ref::MOI.VariableIndex)
    if haskey(m.variable_mapping.refmap, ref.value)
        column = m.variable_mapping.refmap[ref.value]
        if column > 0 && column <= MOI.get(m, MOI.NumberOfVariables())
            return true
        end
    end
    return false
end

#=
    Delete a variable
=#

function MOI.delete!(m::LinQuadSolverInstance, ref::MOI.VariableIndex)
    col = m.variable_mapping[ref]
    lqs_delcols!(m, col, col)
    deleteat!(m.variable_references, col)
    deleteat!(m.variable_primal_solution, col)
    deleteat!(m.variable_dual_solution, col)

    deleteref!(m.variable_mapping, col, ref)
    # deleting from a dict without the key does nothing
    deletebyval!(cmap(m).upper_bound, ref)
    deletebyval!(cmap(m).lower_bound, ref)
    deletebyval!(cmap(m).fixed_bound, ref)
    deletebyval!(cmap(m).interval_bound, ref)

end
MOI.candelete(m::LinQuadSolverInstance, ref::MOI.VariableIndex) = true

# temp fix - change storage for bounds TODO
function deletebyval!(dict::Dict{S,T}, in::T) where {S,T}
    for (key, val) in dict
        if val == in
            delete!(dict, key)
        end
    end
end
#=
    MIP starts
=#
function MOI.set!(m::LinQuadSolverInstance, ::MOI.VariablePrimalStart, ref::MOI.VariableIndex, val::Float64)
    lqs_addmipstarts!(m, [getcol(m, ref)], [val])
end
MOI.canset(m::LinQuadSolverInstance, ::MOI.VariablePrimalStart, ::MOI.VariableIndex) = true

function MOI.set!(m::LinQuadSolverInstance, ::MOI.VariablePrimalStart, refs::Vector{MOI.VariableIndex}, vals::Vector{Float64})
    lqs_addmipstarts!(m, getcol.(m, refs), vals)
end
MOI.canset(m::LinQuadSolverInstance, ::MOI.VariablePrimalStart, ::Vector{MOI.VariableIndex}) = true