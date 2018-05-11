mutable struct LinQuadSOS
    kind::Symbol
    weights::Vector{Float64}
    indices::Vector{Int}
end

mutable struct MockLinQuadModel # <: LinQuadOptInterface.LinQuadOptimizer
    
    sense::Symbol
    
    A::Matrix{Float64}
    b::Vector{Float64}
    c::Vector{Float64}

    #Q::SparseVector{Matrix{Float64}}

    range::Vector{Float64} # constraint UB for range

    lb::Vector{Float64}
    ub::Vector{Float64}

    vartype::Vector{Symbol} # :Bin, :Int, :Con
    contype::Vector{Symbol} # :LEQ :GEQ, :EQ, :RANGE 

    sos::LinQuadSOS{LinQuadSOS}

    termination_status::MOI.TerminationStatusCode
    primal_status::MOI.ResultStatusCode
    dual_status::MOI.ResultStatusCode
 
    variable_primal_solution::Vector{Float64}
    variable_dual_solution::Vector{Float64}

    constraint_primal_solution::Vector{Float64}
    constraint_dual_solution::Vector{Float64}
end

const LQOI = LinQuadOptInterface
const MOI  = LQOI.MOI

const SUPPORTED_OBJECTIVES = [
    LQOI.Linear,
    LQOI.Quad
]

const SUPPORTED_CONSTRAINTS = [
    (LQOI.Linear, LQOI.EQ),
    (LQOI.Linear, LQOI.LE),
    (LQOI.Linear, LQOI.GE),
    (Linear, IV),
    (LQOI.Quad, LQOI.EQ),
    (LQOI.Quad, LQOI.LE),
    (LQOI.Quad, LQOI.GE),
    (LQOI.SinVar, LQOI.EQ),
    (LQOI.SinVar, LQOI.LE),
    (LQOI.SinVar, LQOI.GE),
    (LQOI.SinVar, LQOI.IV),
    (LQOI.SinVar, MOI.ZeroOne),
    (LQOI.SinVar, MOI.Integer),
    (LQOI.VecVar, LQOI.SOS1),
    (LQOI.VecVar, LQOI.SOS2),
    (LQOI.VecVar, MOI.Nonnegatives),
    (LQOI.VecVar, MOI.Nonpositives),
    (LQOI.VecVar, MOI.Zeros),
    (LQOI.VecLin, MOI.Nonnegatives),
    (LQOI.VecLin, MOI.Nonpositives),
    (LQOI.VecLin, MOI.Zeros)
]

mutable struct MockLinQuadOptimizer <: LQOI.LinQuadOptimizer
    LQOI.@LinQuadOptimizerBase
    params::Dict{String,Any}
    MockLinQuadOptimizer(::Void) = new()
end

LQOI.LinearQuadraticModel(::Type{MockLinQuadOptimizer},env) = Model(env::Env,"defaultname")

function MockLinQuadOptimizer(;kwargs...)

    m = MockLinQuadOptimizer(nothing)
    m.params = Dict{String,Any}()
    MOI.empty!(m)
    for (name,value) in kwargs
        m.params[string(name)] = value
        setparam!(m.inner, string(name), value)
    end
    return m
end

function MOI.empty!(m::MockLinQuadOptimizer)
    MOI.empty!(m,m.env)
    for (name,value) in m.params
        setparam!(m.inner, name, value)
    end
end

LQOI.supported_constraints(s::MockLinQuadOptimizer) = SUPPORTED_CONSTRAINTS
LQOI.supported_objectives(s::MockLinQuadOptimizer) = SUPPORTED_OBJECTIVES

cintvec(v::Vector) = convert(Vector{Int32}, v)

LQOI.backend_type(m::MockLinQuadOptimizer, ::MOI.EqualTo{Float64})     = :EQ#Cchar('=')
LQOI.backend_type(m::MockLinQuadOptimizer, ::MOI.LessThan{Float64})    = :LEQ#Cchar('<')
LQOI.backend_type(m::MockLinQuadOptimizer, ::MOI.GreaterThan{Float64}) = :GEQ#Cchar('>')
LQOI.backend_type(m::MockLinQuadOptimizer, ::MOI.Zeros)                = :EQ#Cchar('=')
LQOI.backend_type(m::MockLinQuadOptimizer, ::MOI.Nonpositives)         = :LEQ#Cchar('<')
LQOI.backend_type(m::MockLinQuadOptimizer, ::MOI.Nonnegatives)         = :GEQ#Cchar('>')

# TODO - improve single type
function LQOI.change_variable_bounds!(instance::MockLinQuadOptimizer, colvec, valvec, sensevec)
    lb_len = count(x->x==Cchar('L'), sensevec)
    LB_val = Array{Float64}(0)
    sizehint!(LB_val, lb_len)
    LB_col = Array{Cint}(0)
    sizehint!(LB_col, lb_len)

    ub_len = count(x->x==Cchar('U'), sensevec)
    UB_val = Array{Float64}(0)
    sizehint!(UB_val, ub_len)
    UB_col = Array{Cint}(0)
    sizehint!(UB_col, ub_len)

    for i in eachindex(valvec)
        if sensevec[i] == Cchar('L')
            push!(LB_col, colvec[i])
            push!(LB_val, valvec[i])
        elseif sensevec[i] == Cchar('U')
            push!(UB_col, colvec[i])
            push!(UB_val, valvec[i])
        end
    end

    if lb_len > 0
        for i in eachindex(LB_col)
            instance.inner.lb[LB_col[i]] = LB_val[i]
        end
    end

    if ub_len > 0
        for i in eachindex(UB_col)
            instance.inner.ub[UB_col[i]] = UB_val[i]
        end
    end

    nothing
end

function LQOI.get_variable_lowerbound(instance::MockLinQuadOptimizer, col)
    instance.inner.lb[col]
end

function LQOI.get_variable_upperbound(instance::MockLinQuadOptimizer, col)
    instance.inner.ub[col]
end

function LQOI.get_number_linear_constraints(instance::MockLinQuadOptimizer)
    size(instance.inner.A)[1]
end

function LQOI.add_linear_constraints!(instance::MockLinQuadOptimizer, rowvec, colvec, coefvec, sensevec, rhsvec)
    #add_constrs!(instance.inner, rowvec, colvec, coefvec, sensevec, rhsvec)
end

function LQOI.get_rhs(instance::MockLinQuadOptimizer, row)
    instance.inner.b[row]
end

function LQOI.get_linear_constraint(instance::MockLinQuadOptimizer, idx)
    vals = instance.inner.A[idx,:]
    return (1:length(vals))-1, vals
end

# TODO SPLIT THIS ONE
function LQOI.change_coefficient!(instance::MockLinQuadOptimizer, row, col, coef)
    # if row == 0
    #     set_dblattrlist!(instance.inner, "Obj", Cint[col], Float64[coef])
    # elseif col == 0
    #     set_dblattrlist!(instance.inner, "RHS", Cint[row], Float64[coef])
    # else
    #     chg_coeffs!(instance.inner, row, col, coef)
    #     #TODO fix this function in gurobi
    # end
end

function LQOI.delete_linear_constraints!(instance::MockLinQuadOptimizer, rowbeg, rowend)
    if rowbeg > 1 && rowend < length(instance.inner.b)
        instance.inner.A = vcat(instance.inner.A[1:rowbeg-1,:],instance.inner.A[rowend+1:end,:])
        instance.inner.b = vcat(instance.inner.b[1:rowbeg-1],instance.inner.b[rowend+1:end])
        instance.inner.contype = vcat(instance.inner.contype[1:rowbeg-1],instance.inner.contype[rowend+1:end])
        instance.inner.range = vcat(instance.inner.range[1:rowbeg-1],instance.inner.range[rowend+1:end])
    elseif rowbeg > 1
        instance.inner.A = instance.inner.A[1:rowbeg-1,:]
        instance.inner.b = instance.inner.b[1:rowbeg-1]
        instance.inner.contype = instance.inner.contype[1:rowbeg-1]
        instance.inner.range = instance.inner.range[1:rowbeg-1]
    elseif rowend < length(instance.inner.b)
        instance.inner.A = instance.inner.A[rowend+1:end,:]
        instance.inner.b = instance.inner.b[rowend+1:end]
        instance.inner.contype = instance.inner.contype[rowend+1:end]
        instance.inner.range = instance.inner.range[rowend+1:end]
    else
        instance.inner.A = zeros(0,length(instance.inner.c))
        instance.inner.b = Vector{Float64}[]
        instance.inner.contype = Symbol[]
        instance.inner.range = Vector{Float64}[]
    end
end

# TODO fix types
function variabletype(::MockLinQuadOptimizer, typeval)
    if typeval == "B"
        return :Bin
    elseif typeval == "I"
        return :Int
    else
        return :Con
    end
end
function LQOI.change_variable_types!(instance::MockLinQuadOptimizer, colvec, typevec)
    for i in eachindex(colvec)
        instance.inner.vartype[colvec[i]] = variabletype(instance, typevec[i])
    end
end

function sensetype(::MockLinQuadOptimizer, typeval)
    if typeval == "G"
        return :GEQ
    elseif typeval == "E"
        return :EQ
    elseif typeval == "L"
        return :LEQ
    else
        return :RANGE
    end
end
function LQOI.change_linear_constraint_sense!(instance::MockLinQuadOptimizer, rowvec, sensevec)
    for i in eachindex(rowvec)
        instance.inner.vartype[rowvec[i]] = sensetype(instance, sensevec[i])
    end
end

#=
LQOI.add_sos_constraint!(instance::MockLinQuadOptimizer, colvec, valvec, typ) = (add_sos!(instance.inner, typ, colvec, valvec);update_model!(instance.inner))

LQOI.delete_sos!(instance::MockLinQuadOptimizer, idx1, idx2) = (del_sos!(instance.inner, cintvec(collect(idx1:idx2)));update_model!(instance.inner))

# TODO improve getting processes
function LQOI.get_sos_constraint(instance::MockLinQuadOptimizer, idx)
    A, types = get_sos_matrix(instance.inner)
    line = A[idx,:] #sparse vec
    cols = line.nzind
    vals = line.nzval
    typ = types[idx] == Cint(1) ? :SOS1 : :SOS2
    return cols, vals, typ
end

LQOI.get_number_quadratic_constraints(instance::MockLinQuadOptimizer) = num_qconstrs(instance.inner)

#   NOTE:
# LQOI assumes 0.5 x' Q x, but Gurobi requires x' Q x so we multiply V by 0.5
LQOI.add_quadratic_constraint!(instance::MockLinQuadOptimizer, cols,coefs,rhs,sense, I,J,V) = add_qconstr!(instance.inner, cols, coefs, I, J, 0.5 * V, sense, rhs)

# LQOI.change_range_value!(instance::MockLinQuadOptimizer, rows, vals) = chg_rhsrange!(instance.inner, cintvec(rows), -vals)

# function LQOI.set_quadratic_objective!(instance::MockLinQuadOptimizer, I, J, V)
#     delq!(instance.inner)
#     for i in eachindex(V)
#         if I[i] == J[i]
#             V[i] /= 2
#         end
#     end
#     add_qpterms!(instance.inner, I, J, V)
#     return nothing
# end
=#

function LQOI.set_linear_objective!(instance::MockLinQuadOptimizer, colvec, coefvec)
    nvars = num_vars(instance.inner)
    obj = zeros(Float64, nvars)

    for i in eachindex(colvec)
        obj[colvec[i]] = coefvec[i]
    end

    instance.inner.c = obj
    nothing
end

function LQOI.change_objective_sense!(instance::MockLinQuadOptimizer, symbol)
    if symbol == :min
        instance.inner.sense = :minimize
    else
        instance.inner.sense = :maximize
    end
end

function LQOI.get_linear_objective!(instance::MockLinQuadOptimizer, x)
    copy(instance.inner.c)
end

function LQOI.get_objectivesense(instance::MockLinQuadOptimizer)
    s = instance.inner.sense
    if s == :maximize
        return MOI.MaxSense
    else
        return MOI.MinSense
    end
end

LQOI.get_number_variables(instance::MockLinQuadOptimizer) = length(instance.inner.c)

function LQOI.add_variables!(instance::MockLinQuadOptimizer, int)
    push!(instance.inner.ub,Inf)
    push!(instance.inner.lb,-Inf)
    push!(instance.inner.vartype,:Con)
    instance.inner.A = hcat(instance.inner.A,zeros(length(instance.inner.b)))
end

function LQOI.delete_variables!(instance::MockLinQuadOptimizer, colbeg, colend)
    if colbeg > 1 && colend < length(instance.inner.lb)
        instance.inner.A = hcat(instance.inner.A[:,1:colbeg-1],instance.inner.A[:,colend+1:end])
        instance.inner.lb = vcat(instance.inner.lb[1:colbeg-1],instance.inner.lb[colend+1:end])
        instance.inner.ub = vcat(instance.inner.lb[1:colbeg-1],instance.inner.lb[colend+1:end])
        instance.inner.vartype = vcat(instance.inner.vartype[1:colbeg-1],instance.inner.vartype[colend+1:end])
        instance.inner.range = vcat(instance.inner.range[1:colbeg-1],instance.inner.range[colend+1:end])
    elseif colbeg > 1
        instance.inner.A = instance.inner.A[:,1:colbeg-1]
        instance.inner.lb = instance.inner.lb[1:colbeg-1]
        instance.inner.ub = instance.inner.lb[1:colbeg-1]
        instance.inner.vartype = instance.inner.vartype[1:colbeg-1]
        instance.inner.range = instance.inner.range[1:colbeg-1]
    elseif colend < length(instance.inner.lb)
        instance.inner.A = instance.inner.A[:,colend+1:end]
        instance.inner.lb = instance.inner.lb[colend+1:end]
        instance.inner.ub = instance.inner.lb[colend+1:end]
        instance.inner.vartype = instance.inner.vartype[colend+1:end]
        instance.inner.range = instance.inner.range[colend+1:end]
    else
        instance.inner.A = zeros(length(instance.inner.b),0)
        instance.inner.ub = Vector{Float64}[]
        instance.inner.vartype = Symbol[]
        instance.inner.range = Vector{Float64}[]
    end
end

# function LQOI.add_mip_starts!(instance::MockLinQuadOptimizer, colvec::Vector{Int}, valvec::Vector)
#     x = zeros(num_vars(instance.inner))
#     for (col, val) in zip(colvec, valvec)
#         x[col] = val
#     end
#     loadbasis(instance.inner, x)
# end

LQOI.solve_mip_problem!(instance::MockLinQuadOptimizer) = nothing

LQOI.solve_quadratic_problem!(instance::MockLinQuadOptimizer) = nothing

LQOI.solve_linear_problem!(instance::MockLinQuadOptimizer) = nothing

function LQOI.get_termination_status(instance::MockLinQuadOptimizer)
    return MOI.Success
end

function LQOI.get_primal_status(instance::MockLinQuadOptimizer)
    return MOI.FeasiblePoint
end

function LQOI.get_dual_status(instance::MockLinQuadOptimizer)
    return MOI.FeasiblePoint
end

function LQOI.get_variable_primal_solution!(instance::MockLinQuadOptimizer, place)
    for i in eachindex(place)
        place[i] = instance.inner.variable_primal_solution[i]
    end
end
#=
function LQOI.get_linear_primal_solution!(instance::MockLinQuadOptimizer, place)
    get_dblattrarray!(place, instance.inner, "Slack", 1)
    rhs = get_dblattrarray(instance.inner, "RHS", 1, num_constrs(instance.inner))
    for i in eachindex(place)
        place[i] = -place[i]+rhs[i]
    end
    nothing
end

function LQOI.get_quadratic_primal_solution!(instance::MockLinQuadOptimizer, place)
    get_dblattrarray!(place, instance.inner, "QCSlack", 1)
    rhs = get_dblattrarray(instance.inner, "QCRHS", 1, num_qconstrs(instance.inner))
    for i in eachindex(place)
        place[i] = -place[i]+rhs[i]
    end
    nothing
end

LQOI.get_variable_dual_solution!(instance::MockLinQuadOptimizer, place) = get_dblattrarray!(place, instance.inner, "RC", 1)

LQOI.get_linear_dual_solution!(instance::MockLinQuadOptimizer, place) = get_dblattrarray!(place, instance.inner, "Pi", 1)

LQOI.get_quadratic_dual_solution!(instance::MockLinQuadOptimizer, place) = get_dblattrarray!(place, instance.inner, "QCPi", 1)

LQOI.get_objective_value(instance::MockLinQuadOptimizer) = get_objval(instance.inner)

LQOI.get_objective_bound(instance::MockLinQuadOptimizer) = get_objval(instance.inner)

function LQOI.get_relative_mip_gap(instance::MockLinQuadOptimizer)
    L = get_objval(instance.inner)
    U = get_objbound(instance.inner)
    return abs(U-L)/U
end

LQOI.get_iteration_count(instance::MockLinQuadOptimizer)  = get_iter_count(instance.inner)

LQOI.get_barrier_iterations(instance::MockLinQuadOptimizer) = get_barrier_iter_count(instance.inner)

LQOI.get_node_count(instance::MockLinQuadOptimizer) = get_node_count(instance.inner)

LQOI.get_farkas_dual!(instance::MockLinQuadOptimizer, place) = get_dblattrarray!(place, instance.inner, "FarkasDual", 1)

function hasdualray(instance::MockLinQuadOptimizer)
    try
        get_dblattrarray(instance.inner, "FarkasDual", 1, num_constrs(instance.inner))
        return true
    catch
        return false
    end
end

LQOI.get_unbounded_ray!(instance::MockLinQuadOptimizer, place) = get_dblattrarray!(place, instance.inner, "UnbdRay", 1)

function hasprimalray(instance::MockLinQuadOptimizer)
    try
        get_dblattrarray(instance.inner, "UnbdRay", 1, num_vars(instance.inner))
        return true
    catch
        return false
    end
end

MOI.free!(m::MockLinQuadOptimizer) = free_model(m.inner)
=#