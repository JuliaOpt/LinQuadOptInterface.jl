#=
    This file contains all of the functions that a solver needs to implement
    in order to use LQOI.

       min/max: c'x + x'Qx
    subject to:
                # Variable bounds
                abᵢ <=  xᵢ              <= bbᵢ, i=1,2,...Nx
                # Linear Constraints
                llᵢ <= alᵢ'x             <= ulᵢ, i=1,2,...Nl
                # Quadratic Constraints
                lqᵢ <= aqᵢ' x + x'Qqᵢ' x <= uqᵢ, i=1,2,...Nq
                # SOS1, SOS2 constraints
                # Binary Constraints
                xᵢ ∈ {0, 1}
                # Integer Constraints
                xᵢ ∈ Z
=#

"""
    LinQuadModel(M,env)

Initializes a model given a model type M and a env, that might be a nothing for some solvers.
"""
function LinQuadModel end

# LinQuadSolver # Abstract type
function lqs_setparam!(env, name, val) end
function lqs_setlogfile!(env, path) end
function lqs_getprobtype(m::LinQuadOptimizer) end

lqs_supported_constraints(s) = []
lqs_supported_objectives(s) = []

# Constraints

function lqs_chgbds!(m::LinQuadOptimizer, colvec, valvec, sensevec) end

"""
    lqs_getlb(m, col::Int)::Float64

Get the lower bound of the variable in 1-indexed column `col` of the model `m`.
"""
function lqs_getlb(m::LinQuadOptimizer, col) end

"""
    lqs_getub(m, col::Int)::Float64

Get the upper bound of the variable in 1-indexed column `col` of the model `m`.
"""
function lqs_getub(m::LinQuadOptimizer, col) end

"""
    lqs_getnumrows(m)::Int

Get the number of linear constraints in the model `m`.
"""
function lqs_getnumrows(m::LinQuadOptimizer) end

function lqs_addrows!(m::LinQuadOptimizer, rowvec, colvec, coefvec, sensevec, rhsvec) end

"""
    lqs_getrhs(m, row::Int)::Float64

Get the right-hand side of the linear constraint in the 1-indexed row `row` in
the model `m`.
"""
function lqs_getrhs(m::LinQuadOptimizer, row) end

function lqs_getrows(m::LinQuadOptimizer, rowvec) end

"""
    lqs_getcoef(m, row::Int, col::Int)::Float64

Get the linear coefficient of the variable in column `col`, constraint `row`.
"""
function lqs_getcoef(m::LinQuadOptimizer, row, col) end

"""
    lqs_chgcoef(m, row::Int, col::Int, coef::Float64)::Void

Set the linear coefficient of the variable in column `col`, constraint `row` to
`coef`.
"""
function lqs_chgcoef!(m::LinQuadOptimizer, row, col, coef) end

"""
    lqs_delrows!(m, start_row::Int, end_row::Int)::Void

Delete the linear constraints `start_row`, `start_row+1`, ..., `end_row` from
the model `m`.
"""
function lqs_delrows!(m::LinQuadOptimizer, start_row, end_row) end

function lqs_chgctype!(m::LinQuadOptimizer, colvec, typevec) end

function lqs_chgsense!(m::LinQuadOptimizer, rowvec, sensevec) end

function lqs_vartype_map(m::LinQuadOptimizer) end

function lqs_make_problem_type_integer(m::LinQuadOptimizer) end

function lqs_make_problem_type_continuous(m::LinQuadOptimizer) end

"""
    lqs_addsos!(m, cols::Vector{Int}, vals::Vector{Float64}, typ::Symbol)::Void

Add the SOS constraint to the model `m`. `typ` is either `:SOS1` or `:SOS2`.
"""
function lqs_addsos!(m::LinQuadOptimizer, cols, vals, typ) end

"""
    lqs_delrows!(m, start_idx::Int, end_idx::Int)::Void

Delete the SOS constraints `start_idx`, `start_idx+1`, ..., `end_idx` from
the model `m`.
"""
function lqs_delsos!(m::LinQuadOptimizer, start_idx, end_idx) end

# TODO(@joaquim): what is sertype. Why not a function?
"""
    lqs_sertype_map(m)::Dict

Returns a dictionary that maps  `:SOS1` and `:SOS2` to an appropriate type for
the backend.
"""
function lqs_sertype_map(m::LinQuadOptimizer) end

"""
    lqs_getsos(m, idx::Int)::Tuple{Vector{Int}, Vector{Float64}, Symbol}

Get the SOS constraint `idx` from the model `m`. Returns the triplet
    `(cols, vals, typ)`.
"""
function lqs_getsos(m::LinQuadOptimizer, idx) end

"""
    lqs_getnumqconstrs(m)::Int

Get the number of quadratic constraints in the model `m`.
"""
function lqs_getnumqconstrs(m::LinQuadOptimizer) end

function lqs_addqconstr!(m::LinQuadOptimizer, cols,coefs,rhs,sense, I,J,V) end

function lqs_chgrngval!(m::LinQuadOptimizer, rows, vals) end# later

function lqs_ctrtype_map(m::LinQuadOptimizer) end

#Objective

function lqs_copyquad!(m::LinQuadOptimizer, intvec,intvec2, floatvec) end#?

function lqs_chgobj!(m::LinQuadOptimizer, colvec,coefvec) end

"""
    lqs_chgobjsen(m, sense::Symbol)::Void

Change the optimization sense of the model `m` to `sense`. `sense` must be
`:min` or `:max`.
"""
function lqs_chgobjsen!(m::LinQuadOptimizer, sense) end

#TODO(@joaquimg): why is this not in-place?
"""
    lqs_getobj(m)::Vector{Float64}

Change the linear coefficients of the objective.
"""
function lqs_getobj(m::LinQuadOptimizer) end

"""
    lqs_getobjsen(m)::MOI.OptimizationSense

Get the optimization sense of the model `m`.
"""
function lqs_getobjsen(m::LinQuadOptimizer) end

#Solve

function lqs_mipopt!(m::LinQuadOptimizer) end

function lqs_qpopt!(m::LinQuadOptimizer) end

function lqs_lpopt!(m::LinQuadOptimizer) end

function lqs_getstat(m::LinQuadOptimizer) end

function lqs_solninfo(m::LinQuadOptimizer) end # complex

"""
    lqs_getx!(m, x::Vector{Float64})

Get the primal solution for the variables in the model `m`, and
store in `x`. `x`must have one element for each variable.
"""
function lqs_getx!(m::LinQuadOptimizer, x) end

"""
    lqs_getax!(m, x::Vector{Float64})

Given a set of linear constraints `l <= a'x <= b` in the model `m`, get the
constraint primal `a'x` for each constraint, and store in `x`.
`x` must have one element for each linear constraint.
"""
function lqs_getax!(m::LinQuadOptimizer, x) end

"""
    lqs_getqcax!(m, x::Vector{Float64})

Given a set of quadratic constraints `l <= a'x + x'Qx <= b` in the model `m`,
get the constraint primal `a'x + x'Qx` for each constraint, and store in `x`.
`x` must have one element for each quadratic constraint.
"""
function lqs_getqcax!(m::LinQuadOptimizer, x) end

"""
    lqs_getdj!(m, x::Vector{Float64})

Get the dual solution (reduced-costs) for the variables in the model `m`, and
store in `x`. `x`must have one element for each variable.
"""
function lqs_getdj!(m::LinQuadOptimizer, x) end

"""
    lqs_getpi!(m, x::Vector{Float64})

Get the dual solution for the linear constraints in the model `m`, and
store in `x`. `x`must have one element for each linear constraint.
"""
function lqs_getpi!(m::LinQuadOptimizer, x) end

"""
    lqs_getqcpi!(m, x::Vector{Float64})

Get the dual solution for the quadratic constraints in the model `m`, and
store in `x`. `x`must have one element for each quadratic constraint.
"""
function lqs_getqcpi!(m::LinQuadOptimizer, x) end

"""
    lqs_getobjval!(m)

Get the objective value of the solved model `m`.
"""
function lqs_getobjval(m::LinQuadOptimizer) end

# TODO(@joaquimg): what is this?
function lqs_getbestobjval(m::LinQuadOptimizer) end

"""
    lqs_getmiprelgap!(m)

Get the relative MIP gap of the solved model `m`.
"""
function lqs_getmiprelgap(m::LinQuadOptimizer) end

"""
    lqs_getitcnt!(m)

Get the number of simplex iterations performed during the most recent
optimization of the model `m`.
"""
function lqs_getitcnt(m::LinQuadOptimizer) end

"""
    lqs_getbaritcnt!(m)

Get the number of barrier iterations performed during the most recent
optimization of the model `m`.
"""
function lqs_getbaritcnt(m::LinQuadOptimizer) end

"""
    lqs_getnodecnt!(m)

Get the number of branch-and-cut nodes expolored during the most recent
optimization of the model `m`.
"""
function lqs_getnodecnt(m::LinQuadOptimizer) end

"""
    lqs_dualfarkas!(m, x::Vector{Float64})

Get the farkas dual (certificate of primal infeasiblility) for the linear constraints
in the model `m`, and store in `x`. `x`must have one element for each linear
constraint.
"""
function lqs_dualfarkas!(m::LinQuadOptimizer, x) end

"""
    lqs_getray!(m, x::Vector{Float64})

Get the unbounded ray (certificate of dual infeasiblility) for the linear constraints
in the model `m`, and store in `x`. `x`must have one element for each variable.
"""
function lqs_getray!(m::LinQuadOptimizer, x) end

"""
    lqs_terminationstatus(m)

Get the termination status of the model `m`.
"""
function lqs_terminationstatus(m::LinQuadOptimizer) end

"""
    lqs_primalstatus(m)

Get the primal status of the model `m`.
"""
function lqs_primalstatus(m::LinQuadOptimizer) end

"""
    lqs_dualstatus(m)

Get the dual status of the model `m`.
"""
function lqs_dualstatus(m::LinQuadOptimizer) end

# Variables
"""
    lqs_getnumcols(m)::Int

Get the number of variables in the model `m`.
"""
function lqs_getnumcols(m::LinQuadOptimizer) end

"""
    lqs_newcols!(m, n::Int)::Void

Add `n` new variables to the model `m`.
"""
function lqs_newcols!(m::LinQuadOptimizer, n) end

"""
    lqs_delcols!(m, start_col::Int, end_col::Int)::Void

Delete the columns `start_col`, `start_col+1`, ..., `end_col` from the model `m`.
"""
function lqs_delcols!(m::LinQuadOptimizer, start_col, end_col) end

"""
    lqs_addmipstarts!(m, cols::Vector{Int}, x::Vector{Float64})::Void

Add the MIP start `x` for the variables in the columns `cols` of the model `m`.
"""
function lqs_addmipstarts!(m::LinQuadOptimizer, cols, x) end
