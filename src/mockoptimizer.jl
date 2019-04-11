const LQOI = LinQuadOptInterface
const MOI  = LQOI.MOI

mutable struct LinQuadSOS
    kind::Symbol
    weights::Vector{Float64}
    indices::Vector{Int}
end

mutable struct MockLinQuadModel # <: LinQuadOptInterface.LinQuadOptimizer

    sense::Symbol

    Qcon::Vector{Matrix{Float64}}
    A::Matrix{Float64}
    b::Vector{Float64}
    c::Vector{Float64}
    Qobj::Matrix{Float64}
    objective_constant::Float64
    range::Vector{Float64} # constraint UB for range

    lb::Vector{Float64}
    ub::Vector{Float64}

    vartype::Vector{Cchar} # :Bin, :Int, :Con
    contype::Vector{Cchar} # :LEQ :GEQ, :EQ, :RANGE

    sos::Vector{LinQuadSOS}

    termination_status::MOI.TerminationStatusCode

    primal_status::MOI.ResultStatusCode

    dual_status::MOI.ResultStatusCode

    variable_primal_solution::Vector{Float64}
    variable_dual_solution::Vector{Float64}

    constraint_primal_solution::Vector{Float64}
    constraint_dual_solution::Vector{Float64}

    quadratic_primal_solution::Vector{Float64}
    quadratic_dual_solution::Vector{Float64}

    function MockLinQuadModel(env::Nothing, str::String = "defaultname")
        m = new()

        m.sense = :minimize

        m.Qcon = Matrix{Float64}[]
        m.A = zeros(0,0)
        m.b = zeros(0)
        m.c = zeros(0)
        m.Qobj = zeros(0,0)
        m.objective_constant = 0.0

        m.range = zeros(0)

        m.lb = zeros(0)
        m.ub = zeros(0)

        m.vartype = Cchar[]
        m.contype = Cchar[]

        m.sos = LinQuadSOS[]

        m.termination_status = MOI.OPTIMIZE_NOT_CALLED
        m.primal_status = MOI.NO_SOLUTION
        m.dual_status = MOI.NO_SOLUTION

        m.variable_primal_solution = zeros(0)
        m.variable_dual_solution = zeros(0)

        m.constraint_primal_solution = zeros(0)
        m.constraint_dual_solution = zeros(0)

        m.quadratic_primal_solution = zeros(0)
        m.quadratic_dual_solution = zeros(0)

        return m
    end
end

function Base.show(io::IO, model::MockLinQuadModel)
    println(io, "MockLinQuadModel")
    println(io, "    Sense       $(model.sense)")
    println(io, "    Variables   $(length(model.lb))")
    println(io, "    Constraints $(length(model.b))")
    return
end

"""
Solution builder to create solutions for the mock solver.
Should be called inside `MOI.optimize!(model::LinQuadOptimizer)` right before changing signs. This means that these solutions replicate the low-level solver output.
The idea is to run the test with a solver like Xpress, Gurobi, GLPK or CPLEX to obtain the detailed solution.
"""
function print_low_level_solution(model)
    println()
    println("------------------------------------------")
    println("---------------- solution ----------------")
    println("------------------------------------------")
    println()
    println("---------------- statuses ----------------")
    println()
    @show model.termination_status
    @show model.primal_status
    @show model.dual_status
    println()
    println("----------------- primal -----------------")
    println()
    @show model.variable_primal_solution
    @show model.constraint_primal_solution
    @show model.qconstraint_primal_solution
    println()
    println("------------------ dual- -----------------")
    println()
    @show model.variable_dual_solution
    @show model.constraint_dual_solution
    @show model.qconstraint_dual_solution
    println()
    println("------------------------------------------")
    println()
end


num_vars(inner::MockLinQuadModel) = length(inner.c)
num_cons(inner::MockLinQuadModel) = length(inner.b)

const SUPPORTED_OBJECTIVES = [
    MOI.SingleVariable,
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
    (LQOI.SinVar, MOI.Semiinteger{Float64}),
    (LQOI.SinVar, MOI.Semicontinuous{Float64}),
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
    l_rows::Vector{Int}
    q_rows::Vector{Int}

    termination_status_stored::Vector{MOI.TerminationStatusCode}
    primal_status_stored::Vector{MOI.ResultStatusCode}
    dual_status_stored::Vector{MOI.ResultStatusCode}

    variable_primal_solution_stored::Vector{Vector{Float64}}
    variable_dual_solution_stored::Vector{Vector{Float64}}

    constraint_primal_solution_stored::Vector{Vector{Float64}}
    constraint_dual_solution_stored::Vector{Vector{Float64}}

    quadratic_primal_solution_stored::Vector{Vector{Float64}}
    quadratic_dual_solution_stored::Vector{Vector{Float64}}

    MockLinQuadOptimizer(::Nothing) = new()
end

LQOI.LinearQuadraticModel(::Type{MockLinQuadOptimizer},env::Nothing) = MockLinQuadModel(env,"defaultname")

function MockLinQuadOptimizer(;kwargs...)

    instance = MockLinQuadOptimizer(nothing)
    instance.params = Dict{String,Any}()
    instance.l_rows = Int[]
    instance.q_rows = Int[]

    instance.termination_status_stored = MOI.TerminationStatusCode[]
    instance.primal_status_stored = MOI.ResultStatusCode[]
    instance.dual_status_stored = MOI.ResultStatusCode[]

    instance.variable_primal_solution_stored = Vector{Float64}[]
    instance.variable_dual_solution_stored = Vector{Float64}[]

    instance.constraint_primal_solution_stored = Vector{Float64}[]
    instance.constraint_dual_solution_stored = Vector{Float64}[]

    instance.quadratic_primal_solution_stored = Vector{Float64}[]
    instance.quadratic_dual_solution_stored = Vector{Float64}[]

    MOI.empty!(instance)
    for (name, value) in kwargs
        instance.params[string(name)] = value
        setparam!(instance.inner, string(name), value)
    end
    return instance
end

"""
    unload(from::Vector, to, warn = true)

Helper function to remove the first element of a vector and return it.
If the vector is empty data in `default` is returned instead.
Used in `fakesolve`.
"""
function unload(from, default, warn = true)
    if !isempty(from)
        out = from[1]
        popfirst!(from)
        return out
    else
        if warn
            @warn("No data in the input vector, returning default.")
        end
        return default
    end
end

"""
    fakesolve(instance::MockLinQuadOptimizer)

Set solutions upon solve calls.
Data held in `MockLinQuadOptimizer` stored data is passed to the low-level emulator `MockLinQuadModel`.
"""
function fakesolve(instance::MockLinQuadOptimizer)

    instance.inner.termination_status = unload(instance.termination_status_stored, instance.inner.termination_status)
    instance.inner.primal_status = unload(instance.primal_status_stored, instance.inner.primal_status)
    instance.inner.dual_status = unload(instance.dual_status_stored, instance.inner.dual_status)

    instance.inner.variable_primal_solution = unload(instance.variable_primal_solution_stored, instance.inner.variable_primal_solution)
    instance.inner.variable_dual_solution = unload(instance.variable_dual_solution_stored, instance.inner.variable_dual_solution)

    instance.inner.constraint_primal_solution = unload(instance.constraint_primal_solution_stored, instance.inner.constraint_primal_solution)
    instance.inner.constraint_dual_solution = unload(instance.constraint_dual_solution_stored, instance.inner.constraint_dual_solution)

    instance.inner.quadratic_primal_solution = unload(instance.quadratic_primal_solution_stored, instance.inner.quadratic_primal_solution)
    instance.inner.quadratic_dual_solution = unload(instance.quadratic_dual_solution_stored, instance.inner.quadratic_dual_solution)

    return
end

function set_variable_primal_solution!(instance::MockLinQuadOptimizer, input)
    push!(instance.variable_primal_solution_stored, input)
end
function set_variable_dual_solution!(instance::MockLinQuadOptimizer, input)
    push!(instance.variable_dual_solution_stored, input)
end
function set_constraint_primal_solution!(instance::MockLinQuadOptimizer, input)
    push!(instance.constraint_primal_solution_stored, input)
end
function set_constraint_dual_solution!(instance::MockLinQuadOptimizer, input)
    push!(instance.constraint_dual_solution_stored, input)
end
function set_quadratic_dual_solution!(instance::MockLinQuadOptimizer, input)
    push!(instance.quadratic_dual_solution_stored, input)
end
function set_quadratic_primal_solution!(instance::MockLinQuadOptimizer, input)
    push!(instance.quadratic_primal_solution_stored, input)
end
function set_termination_status!(instance::MockLinQuadOptimizer, input)
    push!(instance.termination_status_stored, input)
end
function set_primal_status!(instance::MockLinQuadOptimizer, input)
    push!(instance.primal_status_stored, input)
end
function set_dual_status!(instance::MockLinQuadOptimizer, input)
    push!(instance.dual_status_stored, input)
end

function set_solution!(instance::MockLinQuadOptimizer;
    variable_primal = Float64[],
    variable_dual = Float64[],
    constraint_primal = Float64[],
    constraint_dual = Float64[],
    quadratic_primal = Float64[],
    quadratic_dual = Float64[],
    termination_status = MOI.OPTIMAL,
    primal_status = MOI.FEASIBLE_POINT,
    dual_status = MOI.FEASIBLE_POINT
    )

    set_variable_primal_solution!(instance, variable_primal)
    set_variable_dual_solution!(instance, variable_dual)
    set_constraint_primal_solution!(instance, constraint_primal)
    set_constraint_dual_solution!(instance, constraint_dual)
    set_quadratic_primal_solution!(instance, quadratic_primal)
    set_quadratic_dual_solution!(instance, quadratic_dual)
    set_termination_status!(instance, termination_status)
    set_primal_status!(instance, primal_status)
    set_dual_status!(instance, dual_status)

    return
end

LQOI.get_number_linear_constraints(instance::MockLinQuadOptimizer) = num_cons(instance.inner) - LQOI.get_number_quadratic_constraints(instance)
function LQOI.get_number_quadratic_constraints(instance::MockLinQuadOptimizer)
    c = 0
    for i in 1:num_cons(instance.inner)
        if !isempty(instance.inner.Qcon[i])
            c += 1
        end
    end
    return c
end

function MOI.empty!(instance::MockLinQuadOptimizer)
    MOI.empty!(instance, nothing)
    # do not empty param
    instance.l_rows = Int[]
    instance.q_rows = Int[]
    for (name, value) in instance.params
        setparam!(instance.inner, name, value)
    end
end

LQOI.supported_constraints(s::MockLinQuadOptimizer) = SUPPORTED_CONSTRAINTS
LQOI.supported_objectives(s::MockLinQuadOptimizer) = SUPPORTED_OBJECTIVES

cintvec(v::Vector) = convert(Vector{Int32}, v)

backend_type(m::MockLinQuadOptimizer, ::MOI.GreaterThan{T}) where T = Cchar('G')
backend_type(m::MockLinQuadOptimizer, ::MOI.LessThan{T}) where T    = Cchar('L')
backend_type(m::MockLinQuadOptimizer, ::MOI.EqualTo{T}) where T     = Cchar('E')
# Implemented separately
# backend_type(m::MockLinQuadOptimizer, ::MOI.Interval{T}) where T    = Cchar('R')

backend_type(m::MockLinQuadOptimizer, ::MOI.Zeros)        = Cchar('E')
backend_type(m::MockLinQuadOptimizer, ::MOI.Nonpositives) = Cchar('L')
backend_type(m::MockLinQuadOptimizer, ::MOI.Nonnegatives) = Cchar('G')

backend_type(m::MockLinQuadOptimizer, ::Val{:Continuous}) = Cchar('C')
backend_type(m::MockLinQuadOptimizer, ::Val{:Upperbound}) = Cchar('U')
backend_type(m::MockLinQuadOptimizer, ::Val{:Lowerbound}) = Cchar('L')

function LQOI.change_variable_bounds!(instance::MockLinQuadOptimizer, colvec, valvec, sensevec)

    lb_len = count(x->x==Cchar('L'), sensevec)
    LB_val = Array{Float64}(undef, 0)
    sizehint!(LB_val, lb_len)
    LB_col = Array{Cint}(undef, 0)
    sizehint!(LB_col, lb_len)

    ub_len = count(x->x==Cchar('U'), sensevec)
    UB_val = Array{Float64}(undef, 0)
    sizehint!(UB_val, ub_len)
    UB_col = Array{Cint}(undef, 0)
    sizehint!(UB_col, ub_len)

    if lb_len + ub_len != length(sensevec)
        error("Invalid values for sensevec")
    end

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
end

function LQOI.get_variable_lowerbound(instance::MockLinQuadOptimizer, col)
    instance.inner.lb[col]
end

function LQOI.get_variable_upperbound(instance::MockLinQuadOptimizer, col)
    instance.inner.ub[col]
end

function LQOI.add_linear_constraints!(instance::MockLinQuadOptimizer, A::CSRMatrix{Float64}, sensevec, rhsvec)

    newrows = length(rhsvec)
    rows = num_cons(instance.inner)
    addedrows = collect((rows+1):(rows+newrows))
    # quadind
    append!(instance.l_rows,addedrows)

    rowvec, colvec, coefvec = A.row_pointers, A.columns, A.coefficients

    rows = length(rhsvec)
    cols = size(instance.inner.A)[2]
    push!(rowvec,length(colvec)+1)

    An = Array(SparseMatrixCSC(cols,rows,rowvec,colvec,coefvec)')

    instance.inner.A = vcat(instance.inner.A,An)

    append!(instance.inner.b,rhsvec)
    append!(instance.inner.contype,sensevec)

    append!(instance.inner.range,[Inf for i in 1:rows])
    append!(instance.inner.Qcon,[zeros(0,0) for i in 1:rows])

    return
end

function LQOI.add_ranged_constraints!(instance::MockLinQuadOptimizer, A::CSRMatrix{Float64}, lowerbound, upperbound)

    newrows = length(lowerbound)
    rows = num_cons(instance.inner)
    addedrows = collect((rows+1):(rows+newrows))
    # quadind
    append!(instance.l_rows,addedrows)

    rowvec, colvec, coefvec = A.row_pointers, A.columns, A.coefficients
    rows = length(lowerbound)
    cols = size(instance.inner.A)[2]
    push!(rowvec,length(colvec)+1)

    # An = Array(sparse(rowvec,colvec,coefvec,rows,cols))
    # @show cols,rows,rowvec,colvec,coefvec
    An = Array(SparseMatrixCSC(cols,rows,rowvec,colvec,coefvec)')

    instance.inner.A = vcat(instance.inner.A,An)

    append!(instance.inner.b,lowerbound)
    append!(instance.inner.contype,['R' for i in 1:rows])

    append!(instance.inner.range,upperbound)
    append!(instance.inner.Qcon,[zeros(0,0) for i in 1:rows])
    return
end

function modify_ranged_constraints!(instance::MockLinQuadOptimizer, rows, lowerbound, upperbound)
    # quadind
    _rows = instance.l_rows[rows]
    for i in _rows
        instance.inner.b[i] = lowerbound[i]
        instance.inner.range[i] = upperbound[i]
    end
end

function LQOI.get_rhs(instance::MockLinQuadOptimizer, row)
    # quadind
    _row = instance.l_rows[row]
    return instance.inner.b[_row]
end
function LQOI.get_quadratic_rhs(instance::MockLinQuadOptimizer, row)
    # quadind
    _row = instance.q_rows[row]
    return instance.inner.b[_row]
end

function LQOI.get_range(instance::MockLinQuadOptimizer, row)
    # quadind
    _row = instance.l_rows[row]
    return instance.inner.b[_row], instance.inner.range[_row]
end

function LQOI.get_linear_constraint(instance::MockLinQuadOptimizer, row)
    # quadind
    _row = instance.l_rows[row]

    vals = instance.inner.A[_row,:]

    outvals = Float64[]
    outinds = Int[]
    for i in eachindex(vals)
        if abs(vals[i]) > 0.0
            push!(outvals,vals[i])
            push!(outinds,i)
        end
    end
    return outinds, outvals
end

function LQOI.get_quadratic_constraint(instance::MockLinQuadOptimizer, row)
    # quadind
    _row = instance.q_rows[row]

    vals = instance.inner.A[_row,:]

    outvals = Float64[]
    outinds = Int[]
    for i in eachindex(vals)
        if abs(vals[i]) > 0.0
            push!(outvals,vals[i])
            push!(outinds,i)
        end
    end

    Q = instance.inner.Qcon[_row]

    #TODO (@joaquim) fix here
    return outinds, outvals, sparse(0.5*Q)
end

function LQOI.change_matrix_coefficient!(instance::MockLinQuadOptimizer, row, col, coef)
    # quadind
    _row = instance.l_rows[row]
    instance.inner.A[_row,col] = coef
end
function LQOI.change_rhs_coefficient!(instance::MockLinQuadOptimizer, row, coef)
    # quadind
    _row = instance.l_rows[row]
    instance.inner.b[_row] = coef
end
function LQOI.change_objective_coefficient!(instance::MockLinQuadOptimizer, col, coef)
    instance.inner.c[col] = coef
end
function LQOI.delete_linear_constraints!(instance::MockLinQuadOptimizer, rowbeg, rowend)
    _rows = instance.l_rows[collect(rowbeg:rowend)]
    survive = setdiff(collect(1:num_cons(instance.inner)), _rows)
    instance.inner.A = instance.inner.A[survive,:]
    instance.inner.b = instance.inner.b[survive]
    instance.inner.contype = instance.inner.contype[survive]
    instance.inner.range = instance.inner.range[survive]
    instance.inner.Qcon = instance.inner.Qcon[survive]
    for i in rowend:-1:rowbeg
        deleteat!(instance.l_rows,i)
    end
end

function LQOI.delete_quadratic_constraints!(instance::MockLinQuadOptimizer, rowbeg, rowend)
    _rows = instance.q_rows[collect(rowbeg:rowend)]
    survive = setdiff(collect(1:num_cons(instance.inner)), _rows)
    instance.inner.A = instance.inner.A[survive,:]
    instance.inner.b = instance.inner.b[survive]
    instance.inner.contype = instance.inner.contype[survive]
    instance.inner.range = instance.inner.range[survive]
    instance.inner.Qcon = instance.inner.Qcon[survive]
    for i in rowend:-1:rowbeg
        deleteat!(instance.q_rows,i)
    end
end

function LQOI.change_variable_types!(instance::MockLinQuadOptimizer, colvec, typevec)
    for i in eachindex(colvec)
        instance.inner.vartype[colvec[i]] = typevec[i]
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
    _rows = instance.l_rows[rowvec]
    for i in eachindex(_rows)
        instance.inner.vartype[_rows[i]] = sensetype(instance, sensevec[i])
    end
end

function LQOI.add_sos_constraint!(instance::MockLinQuadOptimizer, colvec, valvec, typ)
    sos = LinQuadSOS(typ, valvec, colvec)
    push!(instance.inner.sos,sos)
end
function LQOI.delete_sos!(instance::MockLinQuadOptimizer, idx1, idx2)
    delete_between(instance.inner.sos, idx1, idx2)
end
function LQOI.get_sos_constraint(instance::MockLinQuadOptimizer, idx)
    sos = instance.inner.sos[idx]
    return sos.indices, sos.weights, sos.kind
end

function LQOI.add_quadratic_constraint!(instance::MockLinQuadOptimizer, cols, coefs, rhs, sense, I, J, V)
    @assert length(I) == length(J) == length(V)

    nvars = length(instance.inner.c)
    Q = Array(sparse(I,J,V,nvars,nvars))

    a = zeros(nvars)
    for i in eachindex(cols)
        a[cols[i]] = coefs[i]
    end

    instance.inner.A = vcat(instance.inner.A,a')
    push!(instance.inner.b,rhs)
    push!(instance.inner.contype,sense)
    push!(instance.inner.range,Inf)
    push!(instance.inner.Qcon,Q)

    # scalediagonal!(V, I, J, 0.5)
    # scalediagonal!(V, I, J, 2.0)
    push!(instance.q_rows, num_cons(instance.inner))
end


#=
    Objective
=#

function LQOI.set_quadratic_objective!(instance::MockLinQuadOptimizer, I, J, V)
    @assert length(I) == length(J) == length(V)
    n = num_vars(instance.inner)
    Q = Array(sparse(I,J,V,n,n))
    # scalediagonal!(V, I, J, 0.5)
    instance.inner.Qobj = Q
    # scalediagonal!(V, I, J, 2.0)
    return
end

function LQOI.set_linear_objective!(instance::MockLinQuadOptimizer, colvec, coefvec)
    nvars = num_vars(instance.inner)
    obj = zeros(Float64, nvars)

    for i in eachindex(colvec)
        obj[colvec[i]] = coefvec[i]
    end

    instance.inner.c = obj
    return
end

function LQOI.set_constant_objective!(model::MockLinQuadOptimizer, value)
    model.inner.objective_constant = value
    return
end

function LQOI.get_constant_objective(model::MockLinQuadOptimizer)
    return model.inner.objective_constant
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

function LQOI.get_quadratic_terms_objective(instance::MockLinQuadOptimizer)
    I = Int32[]
    J = Int32[]
    V = Float64[]
    Q = instance.inner.Qobj
    n, m = size(Q)
    for i in 1:n, j in 1:n
        if Q[i,j] != 0
            push!(I,i)
            push!(J,j)
            push!(V,Q[i,j])
        end
    end
    return sparse(Q)
end

function LQOI.get_objectivesense(instance::MockLinQuadOptimizer)
    s = instance.inner.sense
    if s == :maximize
        return MOI.MAX_SENSE
    else
        return MOI.MIN_SENSE
    end
end

#=
    Variables
=#

LQOI.get_number_variables(instance::MockLinQuadOptimizer) = length(instance.inner.c)

function LQOI.add_variables!(instance::MockLinQuadOptimizer, int)
    append!(instance.inner.ub,Inf*ones(int))
    append!(instance.inner.lb,-Inf*ones(int))
    append!(instance.inner.c,0.0*ones(int))
    append!(instance.inner.vartype,['C' for i in 1:int])
    instance.inner.A = hcat(instance.inner.A,zeros(length(instance.inner.b),int))
end

function LQOI.delete_variables!(instance::MockLinQuadOptimizer, colbeg, colend)
    instance.inner.Qobj = delete_colrow(instance.inner.Qobj, colbeg, colend)
    instance.inner.Qcon = map(x->delete_colrow(x, colbeg, colend), instance.inner.Qcon)
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
        instance.inner.c = Vector{Float64}[]
        instance.inner.vartype = Cchar[]
    end
end

function delete_colrow(Q, colbeg, colend)
    n, m = size(Q)
    @assert n == m

    if colbeg > n
        # do nothing
        return Q
    elseif colbeg > 1 && colend < n
        Q = hcat(Q[:,1:colbeg-1], Q[:,colend+1:end])
        Q = vcat(Q[1:colbeg-1,:], Q[colend+1:end,:])
        return Q
    elseif colbeg > 1
        return Q[1:colbeg-1, 1:colbeg-1]
    elseif colend < n
        return Q[colend+1:end, colend+1:end]
    else
        return zeros(0,0)
    end
end

function delete_between(V, colbeg, colend)
    n = length(V)
    if colbeg > n
        # do nothing
        return V
    elseif colbeg > 1 && colend < n
        V = vcat(V[1:colbeg-1], V[colend+1:end])
        return V
    elseif colbeg > 1
        return V[1:colbeg-1]
    elseif colend < n
        return V[colend+1:end]
    else
        return zeros(0)
    end
end

#=
    Solve
=#

LQOI.solve_mip_problem!(instance::MockLinQuadOptimizer) = fakesolve(instance)

LQOI.solve_quadratic_problem!(instance::MockLinQuadOptimizer) = fakesolve(instance)

LQOI.solve_linear_problem!(instance::MockLinQuadOptimizer) = fakesolve(instance)

function LQOI.get_termination_status(instance::MockLinQuadOptimizer)
    return instance.inner.termination_status
end

function LQOI.get_primal_status(instance::MockLinQuadOptimizer)
    return instance.inner.primal_status
end

function LQOI.get_dual_status(instance::MockLinQuadOptimizer)
    return instance.inner.dual_status
end

function LQOI.get_variable_primal_solution!(instance::MockLinQuadOptimizer, place)
    for i in eachindex(place)
        place[i] = instance.inner.variable_primal_solution[i]
    end
end

function LQOI.get_linear_primal_solution!(instance::MockLinQuadOptimizer, place)
    for i in eachindex(place)
        place[i] = instance.inner.constraint_primal_solution[i]
    end
end

function LQOI.get_quadratic_primal_solution!(instance::MockLinQuadOptimizer, place)
    for i in eachindex(place)
        place[i] = instance.inner.quadratic_primal_solution[i]
    end
end

function LQOI.get_variable_dual_solution!(instance::MockLinQuadOptimizer, place)
    for i in eachindex(place)
        place[i] = instance.inner.variable_dual_solution[i]
    end
end

function LQOI.get_linear_dual_solution!(instance::MockLinQuadOptimizer, place)
    for i in eachindex(place)
        place[i] = instance.inner.constraint_dual_solution[i]
    end
end

function LQOI.get_quadratic_dual_solution!(instance::MockLinQuadOptimizer, place)
    for i in eachindex(place)
        place[i] = instance.inner.quadratic_dual_solution[i]
    end
end

function get_objval(inner::MockLinQuadModel)
    x = inner.variable_primal_solution
    obj = 0.0
    obj += inner.c' * x
    if !isempty(inner.Qobj)
        Q = Matrix(Symmetric(inner.Qobj, :U))
        obj += (0.5) * x' * Q * x
    end
    return obj + inner.objective_constant
end

LQOI.get_objective_value(instance::MockLinQuadOptimizer) = get_objval(instance.inner)
LQOI.get_objective_bound(instance::MockLinQuadOptimizer) = get_objval(instance.inner)
LQOI.get_relative_mip_gap(instance::MockLinQuadOptimizer) = 0.0
LQOI.get_iteration_count(instance::MockLinQuadOptimizer) = 1
LQOI.get_barrier_iterations(instance::MockLinQuadOptimizer) = 1
LQOI.get_node_count(instance::MockLinQuadOptimizer) = 1

function LQOI.get_farkas_dual!(instance::MockLinQuadOptimizer, place)
    for i in eachindex(place)
        place[i] = instance.inner.constraint_dual_solution[i]
    end
end

function LQOI.get_farkas_dual_bounds!(instance::MockLinQuadOptimizer, place)
    for i in eachindex(place)
        place[i] = instance.inner.variable_dual_solution[i]
    end
end

function hasdualray(instance::MockLinQuadOptimizer)
    return !(NaN in instance.inner.constraint_dual_solution)
end

function LQOI.get_unbounded_ray!(instance::MockLinQuadOptimizer, place)
    for i in eachindex(place)
        place[i] = instance.inner.variable_primal_solution[i]
    end
end

function hasprimalray(instance::MockLinQuadOptimizer)
    return !(NaN in instance.inner.variable_primal_solution)
end

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
function LQOI.add_mip_starts!(instance::MockLinQuadOptimizer, colvec::Vector{Int}, valvec::Vector)
    x = zeros(num_vars(instance.inner))
    for (col, val) in zip(colvec, valvec)
        x[col] = val
    end
    loadbasis(instance.inner, x)
end
=#
