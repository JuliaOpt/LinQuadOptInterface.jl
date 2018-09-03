function set_knapsacktest_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
    termination_status = MOI.Success,
    primal_status = MOI.FeasiblePoint,
    dual_status = MOI.UnknownResultStatus,

    variable_primal = [1.0, 0.0, 0.0, 1.0, 1.0],
    constraint_primal = [9.0],
    quadratic_primal = Float64[],

    variable_dual = [NaN, NaN, NaN, NaN, NaN],
    constraint_dual = [NaN],
    quadratic_dual = Float64[]
    )

    nothing
end

function set_int1test_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
    termination_status = MOI.Success,
    primal_status = MOI.FeasiblePoint,
    dual_status = MOI.UnknownResultStatus,

    variable_primal = [4.0, 5.0, 1.0],
    constraint_primal = [10.0, 15.0],
    quadratic_primal = Float64[],

    variable_dual = [NaN, NaN, NaN],
    constraint_dual = [NaN, NaN],
    quadratic_dual = Float64[]
    )

    nothing
end

function set_int3test_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
    termination_status = MOI.Success,
    primal_status = MOI.FeasiblePoint,
    dual_status = MOI.UnknownResultStatus,

    variable_primal = [1.0, -0.0, -0.0, -0.0, -0.0, 1.0, -0.0, -0.0, -0.0, -0.0, -0.0],
    constraint_primal = [0.9875],
    quadratic_primal = Float64[],

    variable_dual = [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN],
    constraint_dual = [NaN],
    quadratic_dual = Float64[]
    )

    # SOLVE 2
    LQOI.set_solution!(solver,
    termination_status = MOI.Success,
    primal_status = MOI.FeasiblePoint,
    dual_status = MOI.UnknownResultStatus,

    variable_primal = [1.0, -0.0, -0.0, -0.0, -0.0, 1.0, -0.0, -0.0, -0.0, -0.0, -0.0],
    constraint_primal = [0.9875],
    quadratic_primal = Float64[],

    variable_dual = [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN],
    constraint_dual = [NaN],
    quadratic_dual = Float64[],
    )
    nothing
end