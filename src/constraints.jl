#=
    Helper functions to store constraint mappings
=#
cmap(m::LinQuadOptimizer) = m.constraint_mapping
# dummy method for isvalid
constrdict(m::LinQuadOptimizer, ref) = false

_getrhs(set::LE) = set.upper
_getrhs(set::GE) = set.lower
_getrhs(set::EQ) = set.value

function Base.getindex(m::LinQuadOptimizer, c::CI{F,S}) where F where S
    dict = constrdict(m, c)
    return dict[c]
end

function deleteconstraintname!(m::LinQuadOptimizer, ref)
    if haskey(m.constraint_names, ref)
        name = m.constraint_names[ref]
        delete!(m.constraint_names_rev, name)
        delete!(m.constraint_names, ref)
    end
end

"""
    hasinteger(m::LinQuadOptimizer)::Bool

A helper function to determine if the solver instance `m` has any integer
components (i.e. binary, integer, special ordered sets, semicontinuous, or
semi-integer variables).
"""
function hasinteger(m::LinQuadOptimizer)
    (
        length(cmap(m).integer) +
        length(cmap(m).binary) +
        length(cmap(m).sos1) +
        length(cmap(m).sos2) +
        length(cmap(m).semicontinuous) +
        length(cmap(m).semiinteger)
                ) > 0
end

#=
    MOI.isvalid
=#

function MOI.isvalid(m::LinQuadOptimizer, ref::CI{F,S}) where F where S
    dict = constrdict(m, ref)
    if dict == false
        return false
    end
    if haskey(dict, ref)
        return true
    end
    return false
end

#=
    canaddconstraint
=#

function MOI.canaddconstraint(m::LinQuadOptimizer, f::Type{F}, s::Type{S}) where {F<:MOI.AbstractFunction, S<:MOI.AbstractSet}
    return (f,s) in supported_constraints(m)
end

#=
    Get number of constraints
=#

function MOI.canget(m::LinQuadOptimizer, ::MOI.NumberOfConstraints{F, S}) where F where S
    return (F,S) in supported_constraints(m)
end
function MOI.get(m::LinQuadOptimizer, ::MOI.NumberOfConstraints{F, S}) where F where S
    length(constrdict(m, CI{F,S}(UInt(0))))
end

#=
    Get list of constraint references
=#

function MOI.canget(m::LinQuadOptimizer, ::MOI.ListOfConstraintIndices{F, S}) where F where S
    return (F,S) in supported_constraints(m)
end
function MOI.get(m::LinQuadOptimizer, ::MOI.ListOfConstraintIndices{F, S}) where F where S
    sort(collect(keys(constrdict(m, CI{F,S}(UInt(0))))), by=x->x.value)
end

#=
    Get list of constraint types in model
=#

MOI.canget(m::LinQuadOptimizer, ::MOI.ListOfConstraints) = true
function MOI.get(m::LinQuadOptimizer, ::MOI.ListOfConstraints)
    ret = []
    for (F,S) in supported_constraints(m)
        if MOI.get(m, MOI.NumberOfConstraints{F,S}()) > 0
            push!(ret, (F,S))
        end
    end
    ret
end

#=
    Get constraint names
=#

MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintName, ::Type{<:MOI.ConstraintIndex}) = true
function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintName, c::MOI.ConstraintIndex)
    if haskey(m.constraint_names, c)
        m.constraint_names[c]
    else
        ""
    end
end

#=
    Set constraint names
=#

MOI.canset(m::LinQuadOptimizer, ::MOI.ConstraintName, ::Type{<:MOI.ConstraintIndex}) = true
function MOI.set!(m::LinQuadOptimizer, ::MOI.ConstraintName, ref::MOI.ConstraintIndex, name::String)
    if haskey(m.constraint_names_rev, name)
        if m.constraint_names_rev[name] != ref
            error("Duplicate constraint name: $(name)")
        end
    elseif name != ""
        if haskey(m.constraint_names, ref)
            # we're renaming an existing constraint
            old_name = m.constraint_names[ref]
            delete!(m.constraint_names_rev, old_name)
        end
        m.constraint_names[ref] = name
        m.constraint_names_rev[name] = ref
    end
end

#=
    Get constraint by name
=#

# this covers the non-type-stable get(m, ConstraintIndex) case
MOI.canget(m::LinQuadOptimizer, ::Type{MOI.ConstraintIndex}, name::String) = haskey(m.constraint_names_rev, name)
function MOI.get(m::LinQuadOptimizer, ::Type{<:MOI.ConstraintIndex}, name::String)
    m.constraint_names_rev[name]
end

# this covers the type-stable get(m, ConstraintIndex{F,S}, name)::CI{F,S} case
function MOI.canget(m::LinQuadOptimizer, ::Type{FS}, name::String) where FS <: MOI.ConstraintIndex
    haskey(m.constraint_names_rev, name) && typeof(m.constraint_names_rev[name]) == FS
end
function MOI.get(m::LinQuadOptimizer, ::Type{MOI.ConstraintIndex{F,S}}, name::String) where F where S
    m.constraint_names_rev[name]::MOI.ConstraintIndex{F,S}
end


"""
    CSRMatrix{T}

Matrix given in compressed sparse row (CSR) format.

`CSRMatrix` is analgous to the structure in Julia's `SparseMatrixCSC` but with
the rows and columns flipped. It contains three vectors:
 - `row_pointers` is a vector pointing to the start of each row in
    `columns` and `coefficients`;
 - `columns` is a vector of column numbers; and
 - `coefficients` is a vector of corresponding nonzero values.

The length of `row_pointers` is the number of rows in the matrix.

This struct is not a subtype of `AbstractSparseMatrix` as it is intended to be a
collection of the three vectors as they are required by solvers such as Gurobi.
It is not intended to be used for general computation.
"""
struct CSRMatrix{T}
    row_pointers::Vector{Int}
    columns::Vector{Int}
    coefficients::Vector{T}
    function CSRMatrix{T}(row_pointers, columns, coefficients) where T
        @assert length(columns) == length(coefficients)
        new(row_pointers, columns, coefficients)
    end
end

#=
    Below we add constraints.
=#

include("constraints/singlevariable.jl")
include("constraints/vectorofvariables.jl")
include("constraints/scalaraffine.jl")
include("constraints/vectoraffine.jl")
include("constraints/scalarquadratic.jl")
