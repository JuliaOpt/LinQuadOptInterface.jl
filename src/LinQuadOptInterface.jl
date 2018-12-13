__precompile__()
module LinQuadOptInterface

using Compat
using Compat.SparseArrays, Compat.LinearAlgebra

using MathOptInterface
using MathOptInterface.Utilities

const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities

# functions
const Linear = MOI.ScalarAffineFunction{Float64}
const Quad   = MOI.ScalarQuadraticFunction{Float64}
const SinVar = MOI.SingleVariable
const VecVar = MOI.VectorOfVariables
const VecLin = MOI.VectorAffineFunction{Float64}
# sets
const LE = MOI.LessThan{Float64}
const GE = MOI.GreaterThan{Float64}
const EQ = MOI.EqualTo{Float64}
const IV = MOI.Interval{Float64}
const SOS1 = MOI.SOS1{Float64}
const SOS2 = MOI.SOS2{Float64}
# constraint references
const CI{F,S} = MOI.ConstraintIndex{F,S}
const LCI{S} = CI{Linear,S}
const VLCI{S} = CI{VecLin,S}
const QCI{S} = CI{Quad,S}
const SVCI{S} = CI{SinVar,S}
const VVCI{S} = CI{VecVar,S}
# variable reference
const VarInd = MOI.VariableIndex

const LinSets = Union{LE, GE, EQ, IV}
const VecLinSets = Union{MOI.Nonnegatives, MOI.Nonpositives, MOI.Zeros}

struct ConstraintMapping
    # rows in constraint matrix
    less_than::Dict{LCI{LE}, Int}
    greater_than::Dict{LCI{GE}, Int}
    equal_to::Dict{LCI{EQ}, Int}
    interval::Dict{LCI{IV}, Int}

    # vectors of rows in constraint matrix
    nonnegatives::Dict{VLCI{MOI.Nonnegatives}, Vector{Int}}
    nonpositives::Dict{VLCI{MOI.Nonpositives}, Vector{Int}}
    zeros::Dict{VLCI{MOI.Zeros}, Vector{Int}}

    # rows in quadratic constraint matrix
    q_less_than::Dict{QCI{LE}, Int}
    q_greater_than::Dict{QCI{GE}, Int}
    q_equal_to::Dict{QCI{EQ}, Int}

    # references to variable
    upper_bound::Dict{SVCI{LE}, VarInd}
    lower_bound::Dict{SVCI{GE}, VarInd}
    fixed_bound::Dict{SVCI{EQ}, VarInd}
    interval_bound::Dict{SVCI{IV}, VarInd}

    # vectors of rows in constraint matrix
    vv_nonnegatives::Dict{VVCI{MOI.Nonnegatives}, Vector{Int}}
    vv_nonpositives::Dict{VVCI{MOI.Nonpositives}, Vector{Int}}
    vv_zeros::Dict{VVCI{MOI.Zeros}, Vector{Int}}

    integer::Dict{SVCI{MOI.Integer}, VarInd}
    #=
     for some reason CPLEX doesn't respect bounds on a binary variable, so we
     should store the previous bounds so that if we delete the binary constraint
     we can revert to the old bounds
    =#
    binary::Dict{SVCI{MOI.ZeroOne}, Tuple{VarInd, Float64, Float64}}
    sos1::Dict{VVCI{SOS1}, Int}
    sos2::Dict{VVCI{SOS2}, Int}
    semicontinuous::Dict{SVCI{MOI.Semicontinuous{Float64}}, VarInd}
    semiinteger::Dict{SVCI{MOI.Semiinteger{Float64}}, VarInd}

end
ConstraintMapping() = ConstraintMapping(
    Dict{LCI{LE}, Int}(),
    Dict{LCI{GE}, Int}(),
    Dict{LCI{EQ}, Int}(),
    Dict{LCI{IV}, Int}(),
    Dict{VLCI{MOI.Nonnegatives}, Vector{Int}}(),
    Dict{VLCI{MOI.Nonpositives}, Vector{Int}}(),
    Dict{VLCI{MOI.Zeros}, Vector{Int}}(),
    Dict{QCI{LE}, Int}(),
    Dict{QCI{GE}, Int}(),
    Dict{QCI{EQ}, Int}(),
    Dict{SVCI{LE}, VarInd}(),
    Dict{SVCI{GE}, VarInd}(),
    Dict{SVCI{EQ}, VarInd}(),
    Dict{SVCI{IV}, VarInd}(),
    Dict{VVCI{MOI.Nonnegatives}, Vector{Int}}(),
    Dict{VVCI{MOI.Nonpositives}, Vector{Int}}(),
    Dict{VVCI{MOI.Zeros}, Vector{Int}}(),
    Dict{SVCI{MOI.Integer}, VarInd}(),
    Dict{SVCI{MOI.ZeroOne}, Tuple{VarInd, Float64, Float64}}(),
    Dict{VVCI{SOS1}, Int}(),
    Dict{VVCI{SOS2}, Int}(),
    Dict{SVCI{MOI.Semicontinuous{Float64}}, VarInd}(),
    Dict{SVCI{MOI.Semiinteger{Float64}}, VarInd}()
)

"""
    shift_references_after_delete_affine!(m, row)

This function updates all of the references in `m`
after we have deleted row `row` in the affine constraint matrix.
"""
function shift_references_after_delete_affine!(m, row)
    for scalar_affine in [
            cmap(m).less_than,
            cmap(m).greater_than,
            cmap(m).equal_to,
            cmap(m).interval
        ]
        for (key, val) in scalar_affine
            if val > row
                scalar_affine[key] -= 1
            end
        end
    end

    for vector_affine in [
            cmap(m).nonnegatives,
            cmap(m).nonpositives,
            cmap(m).zeros,
            cmap(m).vv_nonnegatives,
            cmap(m).vv_nonpositives,
            cmap(m).vv_zeros
        ]
        for (key, vals) in vector_affine
            for (i, val) in enumerate(vals)
                if val > row
                    vector_affine[key][i] -= 1
                end
            end
        end
    end
end

"""
    shift_references_after_delete_quadratic!(m, row)

This function updates all of the references in `m`
after we have deleted row `row` in the quadratic constraint matrix.
"""
function shift_references_after_delete_quadratic!(m, row)
    for scalar_quadratic in [
            cmap(m).q_less_than,
            cmap(m).q_greater_than,
            cmap(m).q_equal_to
        ]
        for (key, val) in scalar_quadratic
            if val > row
                scalar_quadratic[key] -= 1
            end
        end
    end
end

function Base.isempty(map::ConstraintMapping)

    ret = true
    ret = ret && isempty(map.less_than)
    ret = ret && isempty(map.greater_than)
    ret = ret && isempty(map.equal_to)
    ret = ret && isempty(map.interval)
    ret = ret && isempty(map.nonnegatives)
    ret = ret && isempty(map.nonpositives)
    ret = ret && isempty(map.zeros)
    ret = ret && isempty(map.q_greater_than)
    ret = ret && isempty(map.q_greater_than)
    ret = ret && isempty(map.q_equal_to)
    ret = ret && isempty(map.upper_bound)
    ret = ret && isempty(map.lower_bound)
    ret = ret && isempty(map.fixed_bound)
    ret = ret && isempty(map.interval_bound)
    ret = ret && isempty(map.vv_nonnegatives)
    ret = ret && isempty(map.vv_nonpositives)
    ret = ret && isempty(map.vv_zeros)
    ret = ret && isempty(map.integer)
    ret = ret && isempty(map.binary)
    ret = ret && isempty(map.sos1)
    ret = ret && isempty(map.sos2)
    ret = ret && isempty(map.semiinteger)
    ret = ret && isempty(map.semicontinuous)

    return ret
end

@enum(ObjectiveType,
      SingleVariableObjective,
      AffineObjective,
      QuadraticObjective)

# Abstract + macro
abstract type LinQuadOptimizer <: MOI.AbstractOptimizer end

function Base.show(io::IO, model::LinQuadOptimizer)
    println(io, "A LinQuadOptInterface model with backend:")
    Base.show(io, model.inner)
    return
end

MOIU.supports_default_copy_to(model::LinQuadOptimizer, copy_names::Bool) = true

@enum(VariableType, Continuous, Binary, Integer, Semiinteger, Semicontinuous)

macro LinQuadOptimizerBase(inner_model_type=Any)
    esc(quote
    inner::$inner_model_type

    name::String

    obj_type::LinQuadOptInterface.ObjectiveType
    single_obj_var::Union{Nothing, MOI.VariableIndex}
    obj_sense::MOI.OptimizationSense

    last_variable_reference::UInt64
    variable_mapping::Dict{MOI.VariableIndex, Int}
    variable_names::Dict{MOI.VariableIndex, String}
    variable_names_rev::Dict{String, Set{MOI.VariableIndex}}
    variable_references::Vector{MOI.VariableIndex}
    variable_type::Dict{MOI.VariableIndex, LinQuadOptInterface.VariableType}

    variable_primal_solution::Vector{Float64}
    variable_dual_solution::Vector{Float64}

    last_constraint_reference::UInt64
    constraint_mapping::LinQuadOptInterface.ConstraintMapping

    constraint_constant::Vector{Float64}
    constraint_primal_solution::Vector{Float64}
    constraint_dual_solution::Vector{Float64}

    qconstraint_primal_solution::Vector{Float64}
    qconstraint_dual_solution::Vector{Float64}

    constraint_names::Dict{LinQuadOptInterface.CI, String}
    constraint_names_rev::Dict{String, Set{LinQuadOptInterface.CI}}

    objective_constant::Float64

    termination_status::MOI.TerminationStatusCode
    primal_status::MOI.ResultStatusCode
    dual_status::MOI.ResultStatusCode
    primal_result_count::Int
    dual_result_count::Int

    solvetime::Float64
    end)
end


function MOI.is_empty(m::LinQuadOptimizer)
    ret = true
    ret = ret && m.name == ""
    ret = ret && m.obj_type == AffineObjective
    ret = ret && isa(m.single_obj_var, Nothing)
    ret = ret && m.obj_sense == MOI.MinSense
    ret = ret && m.last_variable_reference == 0
    ret = ret && isempty(m.variable_mapping)
    ret = ret && isempty(m.variable_names)
    ret = ret && isempty(m.variable_names_rev)
    ret = ret && isempty(m.variable_references)
    ret = ret && isempty(m.variable_type)
    ret = ret && isempty(m.variable_primal_solution)
    ret = ret && isempty(m.variable_dual_solution)
    ret = ret && m.last_constraint_reference == 0
    ret = ret && isempty(m.constraint_mapping)
    ret = ret && isempty(m.constraint_constant)
    ret = ret && isempty(m.constraint_primal_solution)
    ret = ret && isempty(m.constraint_dual_solution)
    ret = ret && isempty(m.qconstraint_primal_solution)
    ret = ret && isempty(m.qconstraint_dual_solution)
    ret = ret && isempty(m.constraint_names)
    ret = ret && isempty(m.constraint_names_rev)
    ret = ret && m.objective_constant == 0.0
    ret = ret && m.termination_status == MOI.OptimizeNotCalled
    ret = ret && m.primal_status == MOI.NoSolution
    ret = ret && m.dual_status == MOI.NoSolution
    ret = ret && m.primal_result_count == 0
    ret = ret && m.dual_result_count == 0
    ret = ret && m.solvetime == 0.0

    return ret
end
function MOI.empty!(m::M, env = nothing) where M<:LinQuadOptimizer
    m.name = ""
    m.inner = LinearQuadraticModel(M,env)

    m.obj_type = AffineObjective
    m.single_obj_var = nothing
    # we assume the default is minimization
    m.obj_sense = MOI.MinSense

    m.last_variable_reference = 0
    m.variable_mapping = Dict{MathOptInterface.VariableIndex, Int}()
    m.variable_names = Dict{MathOptInterface.VariableIndex, String}()
    m.variable_names_rev = Dict{String, Set{MathOptInterface.VariableIndex}}()
    m.variable_references = MathOptInterface.VariableIndex[]
    m.variable_type = Dict{MathOptInterface.VariableIndex, VariableType}()
    m.variable_primal_solution = Float64[]
    m.variable_dual_solution = Float64[]

    m.last_constraint_reference = 0
    m.constraint_mapping = LinQuadOptInterface.ConstraintMapping()

    m.constraint_constant = Float64[]
    m.constraint_primal_solution = Float64[]
    m.constraint_dual_solution = Float64[]

    m.qconstraint_primal_solution = Float64[]
    m.qconstraint_dual_solution = Float64[]

    m.constraint_names = Dict{CI, String}()
    m.constraint_names_rev = Dict{String, Set{CI}}()

    m.objective_constant = 0.0

    m.termination_status = MathOptInterface.OptimizeNotCalled
    m.primal_status = MathOptInterface.NoSolution
    m.dual_status = MathOptInterface.NoSolution
    m.primal_result_count = 0
    m.dual_result_count = 0

    m.solvetime = 0.0

    nothing
end

function MOI.get(m::LinQuadOptimizer, ::MOI.Name)
    m.name
end

MOI.supports(::LinQuadOptimizer, ::MOI.Name) = true
function MOI.set(m::LinQuadOptimizer, ::MOI.Name, name::String)
    m.name = name
end

function MOI.supports_constraint(m::LinQuadOptimizer, ft::Type{F}, st::Type{S}) where F <: MOI.AbstractFunction where S <: MOI.AbstractSet
    (ft,st) in supported_constraints(m)
end
function MOI.supports(m::LinQuadOptimizer, ::MOI.ObjectiveFunction{F}) where F <: MOI.AbstractFunction
    F in supported_objectives(m)
end

# a useful helper function
function deleteref!(dict::Dict, i::Int, ref)
    for (key, val) in dict
        if val > i
            dict[key] -= 1
        end
    end
    delete!(dict, ref)
end

include("variables.jl")
include("constraints.jl")
include("objective.jl")
include("solve.jl")
include("copy.jl")

include("solver_interface.jl")

include("mockoptimizer.jl")

end
