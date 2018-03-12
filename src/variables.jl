#=
    Can do
=#

MOI.canaddvariable(::LinQuadOptimizer) = true

#=
    Helper functions
=#
getcol(m::LinQuadOptimizer, ref::VarInd) = m.variable_mapping[ref]
getcol(m::LinQuadOptimizer, v::SinVar) = getcol(m, v.variable)

#=
    Get number of variables
=#

function MOI.get(m::LinQuadOptimizer, ::MOI.NumberOfVariables)
    return lqs_getnumcols(m)
end
MOI.canget(m::LinQuadOptimizer, ::MOI.NumberOfVariables) = true

#=
    Get variable references
=#

function MOI.get(m::LinQuadOptimizer, ::MOI.ListOfVariableIndices)
    return m.variable_references
end
MOI.canget(m::LinQuadOptimizer, ::MOI.ListOfVariableIndices) = true

#=
    Add a variable to the LinQuadOptimizer
Returns a MathOptInterface VariableIndex. So we need to increment the
variable reference counter in the m.variablemapping, and store the column number
in the dictionary
=#

function MOI.addvariable!(m::LinQuadOptimizer)
    lqs_newcols!(m, 1)
    # assumes we add columns linearly
    m.last_variable_reference += 1
    ref = VarInd(m.last_variable_reference)
    m.variable_mapping[ref] = MOI.get(m, MOI.NumberOfVariables())
    push!(m.variable_references, ref)
    push!(m.variable_primal_solution, NaN)
    push!(m.variable_dual_solution, NaN)
    return ref
end

function MOI.addvariables!(m::LinQuadOptimizer, n::Int)
    previous_vars = MOI.get(m, MOI.NumberOfVariables())
    lqs_newcols!(m, n)
    variable_references = VarInd[]
    sizehint!(variable_references, n)
    for i in 1:n
        # assumes we add columns linearly
        m.last_variable_reference += 1
        ref = VarInd(m.last_variable_reference)
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
MOI.candelete(m::LinQuadOptimizer, ref::VarInd) = true

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
function MOI.set!(m::LinQuadOptimizer, ::MOI.VariablePrimalStart, ref::VarInd, val::Float64)
    lqs_addmipstarts!(m, [getcol(m, ref)], [val])
end
MOI.canset(m::LinQuadOptimizer, ::MOI.VariablePrimalStart, ::VarInd) = true

function MOI.set!(m::LinQuadOptimizer, ::MOI.VariablePrimalStart, refs::Vector{VarInd}, vals::Vector{Float64})
    lqs_addmipstarts!(m, getcol.(m, refs), vals)
end
MOI.canset(m::LinQuadOptimizer, ::MOI.VariablePrimalStart, ::Vector{VarInd}) = true