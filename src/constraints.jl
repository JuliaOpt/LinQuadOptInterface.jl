#=
    Helper functions to store constraint mappings
=#
cmap(model::LinQuadOptimizer) = model.constraint_mapping
# dummy method for isvalid
constrdict(::LinQuadOptimizer, ::CI{F, S}) where {F, S} = Dict{CI{F, S}, Any}()

function Base.getindex(model::LinQuadOptimizer, index::CI)
    return constrdict(model, index)[index]
end

# Internal function: delete the constraint name.
function delete_constraint_name(model::LinQuadOptimizer, index::CI)
    if haskey(model.constraint_names, index)
        name = model.constraint_names[index]
        reverse_map = model.constraint_names_rev[name]
        if index in reverse_map
            pop!(reverse_map, index)
        end
        delete!(model.constraint_names, index)
    end
    return
end

"""
    has_integer(model::LinQuadOptimizer)::Bool

A helper function to determine if `model` has any integer components (i.e.
binary, integer, special ordered sets, semicontinuous, or semi-integer
variables).
"""
function has_integer(model::LinQuadOptimizer)
    constraint_map = cmap(model)
    return length(constraint_map.integer) > 0 ||
        length(constraint_map.binary) > 0 ||
        length(constraint_map.sos1) > 0 ||
        length(constraint_map.sos2) > 0 ||
        length(constraint_map.semicontinuous) > 0 ||
        length(constraint_map.semiinteger) > 0
end

#=
    MOI.is_valid
=#

function MOI.is_valid(model::LinQuadOptimizer, index::CI{F,S}) where F where S
    dict = constrdict(model, index)
    return haskey(dict, index)
end

"""
    __assert_valid__(model::LinQuadOptimizer, index::MOI.Index)

Throw an MOI.InvalidIndex error if `MOI.is_valid(model, index) == false`.
"""
function __assert_valid__(model::LinQuadOptimizer, index::MOI.Index)
    if !MOI.is_valid(model, index)
        throw(MOI.InvalidIndex(index))
    end
end

"""
    __assert_supported_constraint__(model::LinQuadOptimizer, ::Type{F}, ::Type{S})

Throw an `UnsupportedConstraint{F, S}` error if the model cannot add constraints
of type `F`-in-`S`.
"""
function __assert_supported_constraint__(model::LinQuadOptimizer, ::Type{F}, ::Type{S}) where {F<:MOI.AbstractFunction, S<:MOI.AbstractSet}
    if !((F,S) in supported_constraints(model))
        throw(MOI.UnsupportedConstraint{F, S}())
    end
end

#=
    Get number of constraints
=#

function MOI.get(model::LinQuadOptimizer, attribute::MOI.NumberOfConstraints{F, S}) where F where S
    __assert_supported_constraint__(model, F, S)
    return length(constrdict(model, CI{F,S}(UInt(0))))
end

#=
    Get list of constraint references
=#

function MOI.get(model::LinQuadOptimizer, attribute::MOI.ListOfConstraintIndices{F, S}) where F where S
    __assert_supported_constraint__(model, F, S)
    dict = constrdict(model, CI{F,S}(UInt(0)))
    indices = collect(keys(dict))
    return sort(indices, by=x->x.value)
end

#=
    Get list of constraint types in model
=#

function MOI.get(model::LinQuadOptimizer, ::MOI.ListOfConstraints)
    return [(F, S) for (F, S) in supported_constraints(model) if
        MOI.get(model, MOI.NumberOfConstraints{F,S}()) > 0]
end

#=
    Get constraint names
=#

function MOI.get(model::LinQuadOptimizer, ::MOI.ConstraintName, index::CI)
    if haskey(model.constraint_names, index)
        return model.constraint_names[index]
    else
        return ""
    end
end
function MOI.supports(
        ::LinQuadOptimizer, ::MOI.ConstraintName, ::Type{<:MOI.ConstraintIndex})
    return true
end

# The rules for MOI are:
# - Allow setting duplicate variable names
# - Throw an error on get if the name is a duplicate
# So we need to store two things.
# 1. a mapping from ConstraintIndex -> Name
# 2. a reverse mapping from Name -> set of ConstraintIndex's with that name
function MOI.set(model::LinQuadOptimizer, ::MOI.ConstraintName,
                  index::MOI.ConstraintIndex, name::String)
    if haskey(model.constraint_names, index)
        # This constraint already has a name, we must be changing it.
        current_name = model.constraint_names[index]
        # Remove `index` from the set of current name.
        pop!(model.constraint_names_rev[current_name], index)
    end
    if name != ""
        # We're changing the name to something non-default, so store it.
        model.constraint_names[index] = name
        if !haskey(model.constraint_names_rev, name)
            model.constraint_names_rev[name] = Set{MOI.ConstraintIndex}()
        end
        push!(model.constraint_names_rev[name], index)
    else
        # We're changing the name to the default, so we don't store it. Note
        # that if `model.constraint_names` doesn't have a key `index`, this does
        # nothing.
        delete!(model.constraint_names, index)
    end
    return
end

#=
    Get constraint by name
=#

function MOI.get(
        model::LinQuadOptimizer, ::Type{<:MOI.ConstraintIndex}, name::String)
    if haskey(model.constraint_names_rev, name)
        if length(model.constraint_names_rev[name]) == 1
            return first(model.constraint_names_rev[name])
        elseif length(model.constraint_names_rev[name]) > 1
            error("Cannot get constraint because the name $(name) is a duplicate.")
        end
    end
    return nothing
end

# this covers the type-stable get(m, ConstraintIndex{F,S}, name)::CI{F,S} case
function MOI.get(model::LinQuadOptimizer, ::Type{MOI.ConstraintIndex{F,S}},
                 name::String) where F where S
    if haskey(model.constraint_names_rev, name)
        if length(model.constraint_names_rev[name]) == 1
            index = first(model.constraint_names_rev[name])
            if isa(index, MOI.ConstraintIndex{F, S})
                return index::MOI.ConstraintIndex{F, S}
            end
        elseif length(model.constraint_names_rev[name]) > 1
            error("Cannot get constraint because the name $(name) is a duplicate.")
        end
    end
    return nothing
end


"""
    CSRMatrix{T}

Matrix given in compressed sparse row (CSR) format.

`CSRMatrix` is analgous to the structure in Julia's `SparseMatrixCSC` but with
the rows and columns flipped. It contains three vectors:
 - `row_pointers` is a vector pointing to the start of each row in
    `columns` and `coefficients` this should end with a `length(coefficients)+1`;
 - `columns` is a vector of column numbers; and
 - `coefficients` is a vector of corresponding nonzero values.

The length of `row_pointers` is the number of rows in the matrix + 1.

This struct is not a subtype of `AbstractSparseMatrix` as it is intended to be a
collection of the three vectors as they are required by some solvers.
Solvers such as Gurobi may need to trim the last index off of row_pointers.
It is not intended to be used for general computation.
"""
struct CSRMatrix{T}
    # This whole definition should really be replaced by and Adjoint{SparseMatrixCSC}
    row_pointers::Vector{Int}
    columns::Vector{Int}
    coefficients::Vector{T}
    function CSRMatrix{T}(row_pointers, columns, coefficients) where T
        @assert length(columns) == length(coefficients)
        @assert length(columns) + 1 == row_pointers[end]
        new(row_pointers, columns, coefficients)
    end
end

#   Move access to an interface so that if we want to change type definition we can.
colvals(mat::CSRMatrix) = mat.columns
row_pointers(mat::CSRMatrix) = mat.row_pointers
row_nonzeros(mat::CSRMatrix) = mat.coefficients

import SparseArrays: sparse
sparse(m::CSRMatrix{T}) where T =sparse(
    SparseMatrixCSC{T,Int}(maximum(m.columns), length(m.row_pointers)-1, m.row_pointers, m.columns, m.coefficients)'
    )


# Lets also depreciate that use of direct property access


#=
    Below we add constraints.
=#

include("constraints/singlevariable.jl")
include("constraints/vectorofvariables.jl")
include("constraints/scalaraffine.jl")
include("constraints/vectoraffine.jl")
include("constraints/scalarquadratic.jl")
