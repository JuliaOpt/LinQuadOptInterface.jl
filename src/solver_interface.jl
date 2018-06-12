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
    backend_type(m, ::MOI.AbstractSet)::CChar

An overloadable type for dispatching the appropriate types to the backends.

For example, GLPK.jl uses `'E'` for `a'x=b` constraints, where as Gurobi.jl uses
`==`.

Three are special cases:
 - `Val{:Continuous}`: for the type of a continuous variable
 - `Val{:Upperbound}`: for the upper bound of a variable
 - `Val{:Lowerbound}`: for the lower bound of a variable

### Defaults

    MOI.GreaterThan  - 'G'
    MOI.LessThan     - 'L'
    MOI.EqualTo      - 'E'

    MOI.Zeros        - 'E'
    MOI.Nonpositives - 'L'
    MOI.Nonnegatives - 'G'

    MOI.ZeroOne        - 'B'
    MOI.Integer        - 'I'
    MOI.Semicontinuous - 'S'
    MOI.Semiinteger    - 'N'

    MOI.SOS1 - :SOS1  # '1'
    MOI.SOS2 - :SOS2  # '2'

    Val{:Continuous} - 'C'
    Val{:Upperbound} - 'U'
    Val{:Lowerbound} - 'L'
"""
backend_type(m::LinQuadOptimizer, ::MOI.GreaterThan{T}) where T = Cchar('G')
backend_type(m::LinQuadOptimizer, ::MOI.LessThan{T}) where T    = Cchar('L')
backend_type(m::LinQuadOptimizer, ::MOI.EqualTo{T}) where T     = Cchar('E')
# Implemented separately
# backend_type(m::LinQuadOptimizer, ::MOI.Interval{T}) where T    = Cchar('R')

backend_type(m::LinQuadOptimizer, ::MOI.Zeros)        = Cchar('E')
backend_type(m::LinQuadOptimizer, ::MOI.Nonpositives) = Cchar('L')
backend_type(m::LinQuadOptimizer, ::MOI.Nonnegatives) = Cchar('G')

backend_type(m::LinQuadOptimizer, ::MOI.ZeroOne) = Cchar('B')
backend_type(m::LinQuadOptimizer, ::MOI.Integer) = Cchar('I')

backend_type(m::LinQuadOptimizer, ::MOI.SOS1{T}) where T = :SOS1  # Cchar('1')
backend_type(m::LinQuadOptimizer, ::MOI.SOS2{T}) where T = :SOS2  # Cchar('2')

backend_type(m::LinQuadOptimizer, ::MOI.Semicontinuous{T}) where T = Cchar('S')
backend_type(m::LinQuadOptimizer, ::MOI.Semiinteger{T}) where T    = Cchar('N')

backend_type(m::LinQuadOptimizer, ::Val{:Continuous}) = Cchar('C')
backend_type(m::LinQuadOptimizer, ::Val{:Upperbound}) = Cchar('U')
backend_type(m::LinQuadOptimizer, ::Val{:Lowerbound}) = Cchar('L')

"""
    LinearQuadraticModel(M::Type{<:LinQuadOptimizer}, env)

Initializes a model given a model type `M` and an `env` that might be a `nothing`
for some solvers.
"""
function LinearQuadraticModel end
@deprecate LinQuadModel LinearQuadraticModel

"""
    supported_constraints(m)::Vector{
        Tuple{MOI.AbstractFunction, MOI.AbstractSet}
    }

Get a list of supported constraint types in the model `m`.

For example, `[(LQOI.Linear, LQOI.EQ)]`
"""
function supported_constraints end
@deprecate lqs_supported_constraints supported_constraints

"""
    supported_objectives(m)::Vector{MOI.AbstractScalarFunction}

Get a list of supported objective types in the model `m`.

For example, `[LQOI.Linear, LQOI.Quad]`
"""
function supported_objectives end
@deprecate lqs_supported_objectives supported_objectives

# Constraints

"""
    change_variable_bounds!(m, cols::Vector{Int}, values::Vector{Float64}, senses::Vector)

Change the bounds of the variable. The sense of the upperbound
is given by `backend_type(m, Val{:Upperbound}())`. The sense
of the lowerbound is given by `backend_type(m, Val{:Lowerbound}())`
"""
function change_variable_bounds! end
@deprecate lqs_chgbds! change_variable_bounds!

"""
    get_variable_lowerbound(m, col::Int)::Float64

Get the lower bound of the variable in 1-indexed column `col` of the model `m`.
"""
function get_variable_lowerbound end
@deprecate lqs_getlb get_variable_lowerbound

"""
    get_variable_upperbound(m, col::Int)::Float64

Get the upper bound of the variable in 1-indexed column `col` of the model `m`.
"""
function get_variable_upperbound end
@deprecate lqs_getub get_variable_upperbound

"""
    get_number_linear_constraints(m)::Int

Get the number of linear constraints in the model `m`.
"""
function get_number_linear_constraints end
@deprecate lqs_getnumrows get_number_linear_constraints

"""
    get_last_linear_constraint_index(m)::Int

Get the index of the last linear constraints added to the model `m`.
"""
function get_last_linear_constraint_index end
get_last_linear_constraint_index(m) = get_number_linear_constraints(m)

"""
    add_linear_constraints!(m, A::CSRMatrix{Float64},
        sense::Vector{Cchar}, rhs::Vector{Float64})::Void

Adds linear constraints of the form `Ax (sense) rhs` to the model `m`.

`sense` and `rhs` contain one element for each row in `A`.
The `sense` is given by `backend_type(m, set)`.

Ranged constraints (`set=MOI.Interval`) should be added via `add_ranged_constraint!`
instead.

See also: `LinQuadOptInterface.CSRMatrix`.
"""
function add_linear_constraints! end
@deprecate lqs_addrows! add_linear_constraints!

"""
    add_ranged_constraints!(m, A::CSRMatrix{Float64},
        lowerbound::Vector{Float64}, upperbound::Vector{Float64})

Adds linear constraints of the form `lowerbound <= Ax <= upperbound` to the
model `m`.

This is a special case compared to standard `add_linear_constraints!` since it
is often implemented via multiple API calls.
"""
function add_ranged_constraints! end

"""
    modify_ranged_constraints!(m, rows::Vector{Int}, lowerbound::Vector{Float64}, upperbound::Vector{Float64})

Modify the lower and upperbounds of a ranged constraint in the model `m`.

This is a special case compared to standard the `change_rhs_coefficient!` since it
is often implemented via multiple API calls.
"""
function modify_ranged_constraints! end

"""
    get_rhs(m, row::Int)::Float64

Get the right-hand side of the linear constraint in the 1-indexed row `row` in
the model `m`.
"""
function get_rhs end
@deprecate lqs_getrhs get_rhs

"""
    get_range(m, row::Int)::Tuple{Float64,Float64}

Get the range which the constraint `row` belongs to. The output of the function is the tuple `lowerbound, upperbound` of bounds: `lowerbound <= a'x < = upperbound`
"""
function get_range end

"""
    get_linear_constraint(m, row::Int)::Tuple{Vector{Int}, Vector{Float64}}

Get the linear component of the constraint in the 1-indexed row `row` in
the model `m`. Returns a tuple of `(cols, vals)`.
"""
function get_linear_constraint end
@deprecate lqs_getrows get_linear_constraint

"""
    change_matrix_coefficient!(m, row, col, coef)

Set the linear coefficient of the variable in column `col`, constraint `row` to
`coef`.
"""
function change_matrix_coefficient! end
@deprecate lqs_chgcoef! change_matrix_coefficient!

"""
change_objective_coefficient!(m, col, coef)

Set the linear coefficient of the variable in column `col` to `coef` in the objective function.
"""
function change_objective_coefficient! end

"""
change_rhs_coefficient!(m, row, coef)

Set the rhs of the constraint in row `row` to `coef`.
"""
function change_rhs_coefficient! end

"""
    delete_linear_constraints!(m, start_row::Int, end_row::Int)::Void

Delete the linear constraints `start_row`, `start_row+1`, ..., `end_row` from
the model `m`.
"""
function delete_linear_constraints! end
@deprecate lqs_delrows! delete_linear_constraints!

"""
    delete_quadratic_constraints!(m, start_row::Int, end_row::Int)::Void

Delete the quadratic constraints `start_row`, `start_row+1`, ..., `end_row` from
the model `m`.
"""
function delete_quadratic_constraints! end

"""
    change_variable_types(m, cols::Vector{Int}, types):Void

Change the variable types. Type is the output of one of:
 - `backend_type(m, ::ZeroOne)`, for binary variables;
 - `backend_type(m, ::Integer)`, for integer variables; and
 - `backend_type(m, Val{:Continuous}())`, for continuous variables.
"""
function change_variable_types! end
@deprecate lqs_chgctype! change_variable_types!

"""
    change_linear_constraint_sense!(m, rows::Vector{Int}, sense::Vector{Symbol})::Void

Change the sense of the linear constraints in `rows` to `sense`.

`sense` is the output of `backend_type(m, set)`, where `set`
is the corresponding set for the row `rows[i]`.

`Interval` constraints require a call to `change_range_value!`.
"""
function change_linear_constraint_sense! end
@deprecate lqs_chgsense! change_linear_constraint_sense!

"""
    make_problem_type_integer(m)::Void

If an explicit call is needed to change the problem type integer (e.g., CPLEX).
"""
function make_problem_type_integer(m::LinQuadOptimizer)
    nothing  # default
end
@deprecate lqs_make_problem_type_integer make_problem_type_integer

"""
    make_problem_type_continuous(m)::Void

If an explicit call is needed to change the problem type continuous (e.g., CPLEX).
"""
function make_problem_type_continuous(m::LinQuadOptimizer)
    nothing  # default
end
@deprecate lqs_make_problem_type_continuous make_problem_type_continuous

"""
    add_sos_constraint!(m, cols::Vector{Int}, vals::Vector{Float64}, typ::Symbol)::Void

Add the SOS constraint to the model `m`. `typ` is either `:SOS1` or `:SOS2`.
"""
function add_sos_constraint! end
@deprecate lqs_addsos! add_sos_constraint!

"""
    delete_sos!(m, start_idx::Int, end_idx::Int)::Void

Delete the SOS constraints `start_idx`, `start_idx+1`, ..., `end_idx` from
the model `m`.
"""
function delete_sos! end
@deprecate lqs_delsos! delete_sos!

"""
    get_sos_constraint(m, idx::Int)::Tuple{Vector{Int}, Vector{Float64}, Symbol}

Get the SOS constraint `idx` from the model `m`. Returns the triplet
    `(cols, vals, typ)`.
"""
function get_sos_constraint end

"""
    get_number_quadratic_constraints(m)::Int

Get the number of quadratic constraints in the model `m`.
"""
function get_number_quadratic_constraints end
@deprecate lqs_getnumqcosntrs get_number_quadratic_constraints

"""
    get_last_quadratic_constraint_index(m)::Int

Get the index of the last quadratic constraints added to the model `m`.
"""
function get_last_quadratic_constraint_index end
get_last_quadratic_constraint_index(m) = get_number_quadratic_constraints(m)

"""
    add_quadratic_constraint!(m, cols::Vector{Int}, coefs::Vector{Float64}, rhs::Float64,
        sense, I::Vector{Int}, J::Vector{Int}, V::Vector{Float64})::Void

Add a quadratic constraint `a'x + 0.5 x' Q x`.
See `add_linear_constraints!` for information of linear component.
Arguments `(I,J,V)` given in triplet form for the Q matrix in `0.5 x' Q x`.
"""
function add_quadratic_constraint! end

"""
    set_quadratic_objective!(m, I::Vector{Int}, J::Vector{Int}, V::Vector{Float64})::Void

Set the quadratic component of the objective. Arguments given in triplet form
for the Q matrix in `0.5 x' Q x`.
"""
function set_quadratic_objective! end
@deprecate lqs_copyquad! set_quadratic_objective!

"""
    get_quadratic_constraint(m, row::Int)::Tuple{Vector{Int}, Vector{Float64}, SparseMatrixCSC{Float64,Int64}}

Get the linear and quadratic components of the constraint in the 1-indexed row `row` in
the model `m`. Returns a tuple of `(lin_cols, lin_vals, Q)`.
Where `Q` represents the matrix in CSC format.
"""
function get_quadratic_constraint end

"""
    get_quadratic_constraint(m, row::Int)::Float64

Get the right hand-side term of quadratic constraint in row `row` in model `m`.
"""
function get_quadratic_rhs end

"""
    set_linear_objective!(m, cols::Vector{Int}, coefs::Vector{Float64})::Void

Set the linear component of the objective.
"""
function set_linear_objective! end
@deprecate lqs_chgobj! set_linear_objective!

"""
    change_objective_sense!(m, sense::Symbol)::Void

Change the optimization sense of the model `m` to `sense`. `sense` must be
`:min` or `:max`.
"""
function change_objective_sense! end
@deprecate lqs_chgobjsen! change_objective_sense!

"""
    get_linear_objective!(m, x::Vector{Float64})

Get the linear coefficients of the objective and store
in `x`.
"""
function get_linear_objective! end
@deprecate lqs_getobj get_linear_objective!

"""
    get_quadratic_terms_objective(m)::SparseMatrixCSC{Float64,Int64}

Get quadratic terms of the objective function returned in sparse CSC format.
"""
function get_quadratic_terms_objective end

"""
    get_objectivesense(m)::MOI.OptimizationSense

Get the optimization sense of the model `m`.
"""
function get_objectivesense end
@deprecate lqs_getobjsen get_objectivesense

"""
    solve_mip_problem!(m)::Void

Solve a mixed-integer model `m`.
"""
function solve_mip_problem! end
@deprecate lqs_mipopt! solve_mip_problem!

"""
    solve_quadratic_problem!(m)::Void

Solve a model `m` with quadratic components.
"""
function solve_quadratic_problem! end
@deprecate lqs_qpopt! solve_quadratic_problem!

"""
    solve_linear_problem!(m)::Void

Solve a linear program `m`.
"""
function solve_linear_problem! end
@deprecate lqs_lpopt! solve_linear_problem!

"""
    get_variable_primal_solution!(m, x::Vector{Float64})

Get the primal solution for the variables in the model `m`, and
store in `x`. `x`must have one element for each variable.
"""
function get_variable_primal_solution! end
@deprecate lqs_getx! get_variable_primal_solution!

"""
    get_linear_primal_solution!(m, x::Vector{Float64})

Given a set of linear constraints `l <= a'x <= b` in the model `m`, get the
constraint primal `a'x` for each constraint, and store in `x`.
`x` must have one element for each linear constraint.
"""
function get_linear_primal_solution! end
@deprecate lqs_getax! get_linear_primal_solution!

"""
    get_quadratic_primal_solution!(m, x::Vector{Float64})

Given a set of quadratic constraints `l <= a'x + x'Qx <= b` in the model `m`,
get the constraint primal `a'x + x'Qx` for each constraint, and store in `x`.
`x` must have one element for each quadratic constraint.
"""
function get_quadratic_primal_solution! end
@deprecate lqs_getqcax! get_quadratic_primal_solution!

"""
    get_variable_dual_solution!(m, x::Vector{Float64})

Get the dual solution (reduced-costs) for the variables in the model `m`, and
store in `x`. `x`must have one element for each variable.
"""
function get_variable_dual_solution! end
@deprecate lqs_getdj! get_variable_dual_solution!

"""
    get_linear_dual_solution!(m, x::Vector{Float64})

Get the dual solution for the linear constraints in the model `m`, and
store in `x`. `x`must have one element for each linear constraint.
"""
function get_linear_dual_solution! end
@deprecate lqs_getpi! get_linear_dual_solution!

"""
    get_quadratic_dual_solution!(m, x::Vector{Float64})

Get the dual solution for the quadratic constraints in the model `m`, and
store in `x`. `x`must have one element for each quadratic constraint.
"""
function get_quadratic_dual_solution! end
@deprecate lqs_getqcpi! get_quadratic_dual_solution!

"""
    get_objective_value(m)

Get the objective value of the solved model `m`.
"""
function get_objective_value end
@deprecate lqs_getobjval get_objective_value

"""
    get_objective_bound(m)

Get the objective bound of the model `m`.
"""
function get_objective_bound end
@deprecate lqs_getbestobjval get_objective_bound

"""
    get_relative_mip_gap(m)

Get the relative MIP gap of the solved model `m`.
"""
function get_relative_mip_gap end
@deprecate lqs_getmiprelgap get_relative_mip_gap

"""
    get_iteration_count(m)

Get the number of simplex iterations performed during the most recent
optimization of the model `m`.
"""
function get_iteration_count end
@deprecate lqs_getitcnt get_iteration_count

"""
    get_barrier_iterations(m)

Get the number of barrier iterations performed during the most recent
optimization of the model `m`.
"""
function get_barrier_iterations end
@deprecate lqs_getbaritcnt get_barrier_iterations

"""
    get_node_count(m)

Get the number of branch-and-cut nodes expolored during the most recent
optimization of the model `m`.
"""
function get_node_count end
@deprecate lqs_getnodecnt get_node_count

"""
    get_farkas_dual!(m, x::Vector{Float64})

Get the farkas dual (certificate of primal infeasiblility) for the linear
constraints in the model `m`, and store in `x`. `x`must have one element for
each linear constraint.
"""
function get_farkas_dual! end
@deprecate lqs_dualfarkas! get_farkas_dual!

"""
    get_unbounded_ray!(m, x::Vector{Float64})

Get the unbounded ray (certificate of dual infeasiblility) for the linear
constraints in the model `m`, and store in `x`. `x`must have one element for
each variable.
"""
function get_unbounded_ray! end
@deprecate lqs_getray! get_unbounded_ray!

"""
    get_termination_status(m)

Get the termination status of the model `m`.
"""
function get_termination_status end
@deprecate lqs_terminationstatus get_termination_status

"""
    get_primal_status(m)

Get the primal status of the model `m`.
"""
function get_primal_status end
@deprecate lqs_primalstatus get_primal_status

"""
    get_dual_status(m)

Get the dual status of the model `m`.
"""
function get_dual_status end
@deprecate lqs_dualstatus get_dual_status

"""
    get_number_variables(m)::Int

Get the number of variables in the model `m`.
"""
function get_number_variables end
@deprecate lqs_getnumcols get_number_variables

"""
    add_variables!(m, n::Int)::Void

Add `n` new variables to the model `m`.
"""
function add_variables! end
@deprecate lqs_newcols! add_variables!

"""
    delete_variables!(m, start_col::Int, end_col::Int)::Void

Delete the columns `start_col`, `start_col+1`, ..., `end_col` from the model `m`.
"""
function delete_variables! end
@deprecate lqs_delcols! delete_variables!

"""
    add_mip_starts!(m, cols::Vector{Int}, x::Vector{Float64})::Void

Add the MIP start `x` for the variables in the columns `cols` of the model `m`.
"""
function add_mip_starts! end
@deprecate lqs_addmipstarts! add_mip_starts!
