function set_knapsacktest_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.UnknownResultStatus)

    LQOI.set_variable_primal_solution!(solver, [1.0, 0.0, 0.0, 1.0, 1.0])
    LQOI.set_constraint_primal_solution!(solver, [9.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [NaN, NaN, NaN, NaN, NaN])
    LQOI.set_constraint_dual_solution!(solver, [NaN])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    nothing
end

function set_int1test_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.UnknownResultStatus)

    LQOI.set_variable_primal_solution!(solver, [4.0, 5.0, 1.0])
    LQOI.set_constraint_primal_solution!(solver, [10.0, 15.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [NaN, NaN, NaN])
    LQOI.set_constraint_dual_solution!(solver, [NaN, NaN])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    nothing
end

function set_int3test_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.UnknownResultStatus)

    LQOI.set_variable_primal_solution!(solver, [1.0, -0.0, -0.0, -0.0, -0.0, 1.0, -0.0, -0.0, -0.0, -0.0, -0.0])
    LQOI.set_constraint_primal_solution!(solver, [0.9875])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN])
    LQOI.set_constraint_dual_solution!(solver, [NaN])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    # SOLVE 2
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.UnknownResultStatus)

    LQOI.set_variable_primal_solution!(solver, [1.0, -0.0, -0.0, -0.0, -0.0, 1.0, -0.0, -0.0, -0.0, -0.0, -0.0])
    LQOI.set_constraint_primal_solution!(solver, [0.9875])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN])
    LQOI.set_constraint_dual_solution!(solver, [NaN])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])
    nothing
end