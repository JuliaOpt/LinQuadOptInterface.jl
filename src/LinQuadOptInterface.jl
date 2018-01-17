module LinQuadOptInterface

using MathOptInterface

const MOI = MathOptInterface

# functions
const Linear = MOI.ScalarAffineFunction{Float64}
const Quad   = MOI.ScalarQuadraticFunction{Float64}
const SinVar = MOI.SingleVariable
const VecVar = MOI.VectorOfVariables
const VecLin = MOI.VectorAffineFunction{Float64}
# sets
const LE     = MOI.LessThan{Float64}
const GE     = MOI.GreaterThan{Float64}
const EQ     = MOI.EqualTo{Float64}
const IV     = MOI.Interval{Float64}
# constraint references
const CI{F,S} = MOI.ConstraintIndex{F,S}
const LCI{S} = CI{Linear,S}
const VLCI{S} = CI{VecLin,S}
const QCI{S} = CI{Quad,S}
const SVCI{S} = CI{SinVar, S}
const VVCI{S} = CI{VecVar, S}
# variable reference
const VarInd = MOI.VariableIndex

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
    interval_bound::Dict{SVCI{MOI.Interval{Float64}}, VarInd}

    # vectors of rows in constraint matrix
    vv_nonnegatives::Dict{VVCI{MOI.Nonnegatives}, Vector{VarInd}}
    vv_nonpositives::Dict{VVCI{MOI.Nonpositives}, Vector{VarInd}}
    vv_zeros::Dict{VVCI{MOI.Zeros}, Vector{VarInd}}

    integer::Dict{SVCI{MOI.Integer}, VarInd}
    #=
     for some reason CPLEX doesn't respect bounds on a binary variable, so we
     should store the previous bounds so that if we delete the binary constraint
     we can revert to the old bounds
    =#
    binary::Dict{SVCI{MOI.ZeroOne}, Tuple{VarInd, Float64, Float64}}
    sos1::Dict{VVCI{MOI.SOS1}, Int}
    sos2::Dict{VVCI{MOI.SOS2}, Int}
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
    Dict{VVCI{MOI.Nonnegatives}, Vector{VarInd}}(),
    Dict{VVCI{MOI.Nonpositives}, Vector{VarInd}}(),
    Dict{VVCI{MOI.Zeros}, Vector{VarInd}}(),
    Dict{SVCI{MOI.Integer}, VarInd}(),
    Dict{SVCI{MOI.ZeroOne}, Tuple{VarInd, Float64, Float64}}(),
    Dict{VVCI{MOI.SOS1}, Int}(),
    Dict{VVCI{MOI.SOS2}, Int}()
)
macro def(name, definition)
    return quote
        macro $(esc(name))()
            esc($(Expr(:quote, definition)))
        end
    end
end

# Abstract + macro
abstract type LinQuadSolverInstance <: MOI.AbstractSolverInstance end
@def LinQuadSolverInstanceBase begin
    inner::Model

    obj_is_quad::Bool

    last_variable_reference::UInt64
    variable_mapping::Dict{MathOptInterface.VariableIndex, Int}
    variable_references::Vector{MathOptInterface.VariableIndex}

    variable_primal_solution::Vector{Float64}
    variable_dual_solution::Vector{Float64}

    last_constraint_reference::UInt64
    constraint_mapping::LinQuadOptInterface.ConstraintMapping

    constraint_constant::Vector{Float64}

    constraint_primal_solution::Vector{Float64}
    constraint_dual_solution::Vector{Float64}

    qconstraint_primal_solution::Vector{Float64}
    qconstraint_dual_solution::Vector{Float64}

    objective_constant::Float64

    termination_status::MathOptInterface.TerminationStatusCode
    primal_status::MathOptInterface.ResultStatusCode
    dual_status::MathOptInterface.ResultStatusCode
    primal_result_count::Int
    dual_result_count::Int

    solvetime::Float64
end

# function MOI.supportsproblem(s::LinQuadSolverInstance, objective_type, constraint_types)
#     if !(objective_type in lqs_supported_objectives(s))
#         return false
#     end
#     for c in constraint_types
#         if !(c in lqs_supported_constraints(s))
#             return false
#         end
#     end
#     return true
# end

@def LinQuadSolverInstanceBaseInit begin
    Model(env),
    false,
    0,
    Dict{MathOptInterface.VariableIndex, Int}(),
    MathOptInterface.VariableIndex[],
    Float64[],
    Float64[],
    0,
    LinQuadOptInterface.ConstraintMapping(),
    Float64[],
    Float64[],
    Float64[],
    Float64[],
    Float64[],
    0.0,
    MathOptInterface.OtherError, # not solved
    MathOptInterface.UnknownResultStatus,
    MathOptInterface.UnknownResultStatus,
    0,
    0,
    0.0
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

# function problemtype(m::LinQuadSolverInstance)
#     code = lqs_getprobtype(m)
#     PROB_TYPE_MAP[code]
# end

include("lqoi_variables.jl")
include("lqoi_constraints.jl")
include("lqoi_objective.jl")
include("lqoi_solve.jl")

include("ref.jl")

end