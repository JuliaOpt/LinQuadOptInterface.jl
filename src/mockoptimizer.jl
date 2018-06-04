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

    vartype::Vector{Cchar} # :Bin, :Int, :Con
    contype::Vector{Cchar} # :LEQ :GEQ, :EQ, :RANGE

    # sos::LinQuadSOS{LinQuadSOS}

    termination_status::MOI.TerminationStatusCode
    primal_status::MOI.ResultStatusCode
    dual_status::MOI.ResultStatusCode

    variable_primal_solution::Vector{Float64}
    variable_dual_solution::Vector{Float64}

    constraint_primal_solution::Vector{Float64}
    constraint_dual_solution::Vector{Float64}
    function MockLinQuadModel(env::Void,str::String="defaultname")
        m = new()

        m.sense = :minimize

        m.A = zeros(0,0)
        m.b = zeros(0)
        m.c = zeros(0)

        m.range = zeros(0)

        m.lb = zeros(0)
        m.ub = zeros(0)

        m.vartype = Cchar[]
        m.contype = Cchar[]

        m.termination_status = MOI.Success
        m.primal_status = MOI.FeasiblePoint
        m.dual_status = MOI.FeasiblePoint

        m.variable_primal_solution = zeros(0)
        m.variable_dual_solution = zeros(0)

        m.constraint_primal_solution = zeros(0)
        m.constraint_dual_solution = zeros(0)

        return m
    end
end
num_vars(inner::MockLinQuadModel) = length(inner.c)
num_cons(inner::MockLinQuadModel) = length(inner.rhs)

function set_variable_primal_solution!(inner::MockLinQuadModel,input::Vector)
    m.variable_primal_solution = input
end
function set_variable_dual_solution!(inner::MockLinQuadModel,input::Vector)
    m.variable_dual_solution = input
end
function set_constraint_primal_solution!(inner::MockLinQuadModel,input::Vector)
    m.constraint_primal_solution = input
end
function set_constraint_dual_solution!(inner::MockLinQuadModel,input::Vector)
    m.constraint_dual_solution = input
end
function set_termination_status!(inner::MockLinQuadModel,input)
    m.termination_status = input
end
function set_primal_status!(inner::MockLinQuadModel,input)
    m.primal_status = input
end
function set_dual_status!(inner::MockLinQuadModel,input)
    m.dual_status = input
end

const LQOI = LinQuadOptInterface
const MOI  = LQOI.MOI

const SUPPORTED_OBJECTIVES = [
    LQOI.Linear,
    #LQOI.Quad
]

const SUPPORTED_CONSTRAINTS = [
    (LQOI.Linear, LQOI.EQ),
    (LQOI.Linear, LQOI.LE),
    (LQOI.Linear, LQOI.GE),
    (Linear, IV),
    #(LQOI.Quad, LQOI.EQ),
    #(LQOI.Quad, LQOI.LE),
    #(LQOI.Quad, LQOI.GE),
    (LQOI.SinVar, LQOI.EQ),
    (LQOI.SinVar, LQOI.LE),
    (LQOI.SinVar, LQOI.GE),
    (LQOI.SinVar, LQOI.IV),
    (LQOI.SinVar, MOI.ZeroOne),
    (LQOI.SinVar, MOI.Integer),
    #(LQOI.VecVar, LQOI.SOS1),
    #(LQOI.VecVar, LQOI.SOS2),
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

LQOI.LinearQuadraticModel(::Type{MockLinQuadOptimizer},env::Void) = MockLinQuadModel(env,"defaultname")

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
    MOI.empty!(m,nothing)
    for (name,value) in m.params
        setparam!(m.inner, name, value)
    end
end

LQOI.supported_constraints(s::MockLinQuadOptimizer) = SUPPORTED_CONSTRAINTS
LQOI.supported_objectives(s::MockLinQuadOptimizer) = SUPPORTED_OBJECTIVES

cintvec(v::Vector) = convert(Vector{Int32}, v)

LQOI.backend_type(m::MockLinQuadOptimizer, ::MOI.EqualTo{Float64})     = Cchar('E')#Cchar('=')
LQOI.backend_type(m::MockLinQuadOptimizer, ::MOI.LessThan{Float64})    = Cchar('L')#Cchar('<')
LQOI.backend_type(m::MockLinQuadOptimizer, ::MOI.GreaterThan{Float64}) = Cchar('G')#Cchar('>')
LQOI.backend_type(m::MockLinQuadOptimizer, ::MOI.Zeros)                = Cchar('E')#Cchar('=')
LQOI.backend_type(m::MockLinQuadOptimizer, ::MOI.Nonpositives)         = Cchar('L')#Cchar('<')
LQOI.backend_type(m::MockLinQuadOptimizer, ::MOI.Nonnegatives)         = Cchar('G')#Cchar('>')

# TODO - improve single type
function LQOI.change_variable_bounds!(instance::MockLinQuadOptimizer, colvec, valvec, sensevec)
    # colvec = colvec+1

    lb_len = count(x->x==Cchar('L'), sensevec)
    LB_val = Array{Float64}(0)
    sizehint!(LB_val, lb_len)
    LB_col = Array{Cint}(0)
    sizehint!(LB_col, lb_len)

    ub_len = count(x->x==Cchar('G'), sensevec)
    UB_val = Array{Float64}(0)
    sizehint!(UB_val, ub_len)
    UB_col = Array{Cint}(0)
    sizehint!(UB_col, ub_len)

    for i in eachindex(valvec)
        if sensevec[i] == Cchar('L')
            push!(LB_col, colvec[i])
            push!(LB_val, valvec[i])
        elseif sensevec[i] == Cchar('G')
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

function LQOI.add_linear_constraints!(instance::MockLinQuadOptimizer, A::CSRMatrix{Float64}, sensevec, rhsvec)
    rowvec, colvec, coefvec = A.row_pointers, A.column_indices, A.data

    rows = length(rhsvec)
    cols = size(instance.inner.A)[2]
    push!(rowvec,length(colvec)+1)

    # An = full(sparse(rowvec,colvec,coefvec,rows,cols))
    # @show cols,rows,rowvec,colvec,coefvec
    An = full(SparseMatrixCSC(cols,rows,rowvec,colvec,coefvec)')

    instance.inner.A = vcat(instance.inner.A,An)

    append!(instance.inner.b,rhsvec)
    append!(instance.inner.contype,sensevec)

    append!(instance.inner.range,[Inf for i in 1:rows])
    nothing
end

function LQOI.add_ranged_constraints!(instance::MockLinQuadOptimizer, A::CSRMatrix{Float64}, rhsvec, ubvec)
    rowvec, colvec, coefvec = A.row_pointers, A.column_indices, A.data
    rows = length(rhsvec)
    cols = size(instance.inner.A)[2]
    push!(rowvec,length(colvec)+1)

    # An = full(sparse(rowvec,colvec,coefvec,rows,cols))
    # @show cols,rows,rowvec,colvec,coefvec
    An = full(SparseMatrixCSC(cols,rows,rowvec,colvec,coefvec)')

    instance.inner.A = vcat(instance.inner.A,An)

    append!(instance.inner.b,rhsvec)
    append!(instance.inner.contype,['R' for i in 1:rows])

    append!(instance.inner.range,ubvec)

    nothing
end

function modify_ranged_constraints!(instance::MockLinQuadOptimizer, rows, rhsvec, ubvec)
    for i in rows
        instance.inner.b[i] = rhsvec[i]
        instance.inner.range[i] = ubvec[i]
    end
    nothing
end

function LQOI.get_rhs(instance::MockLinQuadOptimizer, row)
    instance.inner.b[row]
end

function LQOI.get_linear_constraint(instance::MockLinQuadOptimizer, idx)
    vals = instance.inner.A[idx,:]

    outvals = Float64[]
    outinds = Int[]
    for i in eachindex(vals)
        if abs(vals[i]) > 0.0
            push!(outvals,vals[i])
            push!(outinds,i-1)
        end
    end
    return outinds, outvals
end

# TODO SPLIT THIS ONE
function LQOI.change_matrix_coefficient!(instance::MockLinQuadOptimizer, row, col, coef)
    # if row == 0
    #     set_dblattrlist!(instance.inner, "Obj", Cint[col], Float64[coef])
    # elseif col == 0
    #     set_dblattrlist!(instance.inner, "RHS", Cint[row], Float64[coef])
    # else
    #     chg_coeffs!(instance.inner, row, col, coef)
    #     #TODO fix this function in gurobi
    # end
end
function LQOI.change_rhs_coefficient!(instance::MockLinQuadOptimizer, row, coef)
end
function LQOI.change_objective_coefficient!(instance::MockLinQuadOptimizer, col, coef)
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
        instance.inner.contype = Cchar[]
        instance.inner.range = Vector{Float64}[]
    end
end

# TODO fix types
function variabletype(::MockLinQuadOptimizer, typeval)
    if typeval == 'B'
        return 'B'
    elseif typeval == 'I'
        return 'I'
    else
        return 'C'
    end
end
function LQOI.change_variable_types!(instance::MockLinQuadOptimizer, colvec, typevec)
    for i in eachindex(colvec)
        instance.inner.vartype[colvec[i]] = variabletype(instance, typevec[i])
    end
end

function sensetype(::MockLinQuadOptimizer, typeval)
    if typeval == 'G'
        return 'G'
    elseif typeval == 'E'
        return 'E'
    elseif typeval == 'L'
        return 'L'
    else
        return 'R'
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
    for i in eachindex(x)
        x[i] = instance.inner.c[i]
    end
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
    append!(instance.inner.ub,Inf*ones(int))
    append!(instance.inner.lb,-Inf*ones(int))
    append!(instance.inner.c,0.0*ones(int))
    append!(instance.inner.vartype,['C' for i in 1:int])
    instance.inner.A = hcat(instance.inner.A,zeros(length(instance.inner.b),int))
end

function LQOI.delete_variables!(instance::MockLinQuadOptimizer, colbeg, colend)
    if colbeg > 1 && colend < length(instance.inner.lb)
        instance.inner.A = hcat(instance.inner.A[:,1:colbeg-1],instance.inner.A[:,colend+1:end])
        instance.inner.lb = vcat(instance.inner.lb[1:colbeg-1],instance.inner.lb[colend+1:end])
        instance.inner.ub = vcat(instance.inner.ub[1:colbeg-1],instance.inner.ub[colend+1:end])
        instance.inner.c = vcat(instance.inner.c[1:colbeg-1],instance.inner.c[colend+1:end])
        instance.inner.vartype = vcat(instance.inner.vartype[1:colbeg-1],instance.inner.vartype[colend+1:end])
    elseif colbeg > 1
        instance.inner.A = instance.inner.A[:,1:colbeg-1]
        instance.inner.lb = instance.inner.lb[1:colbeg-1]
        instance.inner.ub = instance.inner.ub[1:colbeg-1]
        instance.inner.c = instance.inner.c[1:colbeg-1]
        instance.inner.vartype = instance.inner.vartype[1:colbeg-1]
    elseif colend < length(instance.inner.lb)
        instance.inner.A = instance.inner.A[:,colend+1:end]
        instance.inner.lb = instance.inner.lb[colend+1:end]
        instance.inner.ub = instance.inner.ub[colend+1:end]
        instance.inner.c = instance.inner.c[colend+1:end]
        instance.inner.vartype = instance.inner.vartype[colend+1:end]
    else
        instance.inner.A = zeros(length(instance.inner.b),0)
        instance.inner.lb = Vector{Float64}[]
        instance.inner.ub = Vector{Float64}[]
        instance.inner.cß = Vector{Float64}[]
        instance.inner.vartype = Cchar[]
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
    nothing
end

function LQOI.get_linear_primal_solution!(instance::MockLinQuadOptimizer, place)
    # Fix QUADRATICS
    for i in eachindex(place)
        place[i] = instance.inner.constraint_primal_solution[i]
    end
    nothing
end

#=
function LQOI.get_quadratic_primal_solution!(instance::MockLinQuadOptimizer, place)
    get_dblattrarray!(place, instance.inner, "QCSlack", 1)
    rhs = get_dblattrarray(instance.inner, "QCRHS", 1, num_qconstrs(instance.inner))
    for i in eachindex(place)
        place[i] = -place[i]+rhs[i]
    end
    nothing
end
=#

function LQOI.get_variable_dual_solution!(instance::MockLinQuadOptimizer, place)
    for i in eachindex(place)
        place[i] = instance.inner.variable_dual_solution[i]
    end
    nothing
end

function LQOI.get_linear_dual_solution!(instance::MockLinQuadOptimizer, place)
    # Fix QUADRATICS
    for i in eachindex(place)
        place[i] = instance.inner.constraint_dual_solution[i]
    end
    nothing
end

#LQOI.get_quadratic_dual_solution!(instance::MockLinQuadOptimizer, place) = get_dblattrarray!(place, instance.inner, "QCPi", 1)

LQOI.get_objective_value(instance::MockLinQuadOptimizer) = get_objval(instance.inner)

#=
LQOI.get_objective_bound(instance::MockLinQuadOptimizer) = get_objval(instance.inner)

function LQOI.get_relative_mip_gap(instance::MockLinQuadOptimizer)
    L = get_objval(instance.inner)
    U = get_objbound(instance.inner)
    return abs(U-L)/U
end

LQOI.get_iteration_count(instance::MockLinQuadOptimizer)  = get_iter_count(instance.inner)

LQOI.get_barrier_iterations(instance::MockLinQuadOptimizer) = get_barrier_iter_count(instance.inner)

LQOI.get_node_count(instance::MockLinQuadOptimizer) = get_node_count(instance.inner)
=#
#=
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
=#
MOI.free!(m::MockLinQuadOptimizer) = nothing
