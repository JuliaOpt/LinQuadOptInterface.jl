function set_lin1test_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [1.0, 0.0, 2.0])
    LQOI.set_constraint_primal_solution!(solver, [1.0, 0.0, 2.0, 3.0, 2.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [0.0, 0.0, 0.0])
    LQOI.set_constraint_dual_solution!(solver, [-0.0, 2.0, -0.0, -3.0, -1.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    nothing
end

function set_lin2test_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [-4.0, -3.0, 16.0, 0.0])
    LQOI.set_constraint_primal_solution!(solver, [-4.0, -3.0, 12.0, -3.0, 16.0, 0.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [0.0, 0.0, 0.0, 0.0])
    LQOI.set_constraint_dual_solution!(solver, [7.0, 2.0, -4.0, -0.0, -0.0, 7.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    nothing
end

function set_lin3test_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.UnknownResultStatus)
    LQOI.set_dual_status!(solver, MOI.InfeasibilityCertificate)

    LQOI.set_variable_primal_solution!(solver, [NaN])
    LQOI.set_constraint_primal_solution!(solver, [NaN, NaN])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [NaN])
    LQOI.set_constraint_dual_solution!(solver, [1.0, -1.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    nothing
end

function set_lin3test_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.UnknownResultStatus)
    LQOI.set_dual_status!(solver, MOI.InfeasibilityCertificate)

    LQOI.set_variable_primal_solution!(solver, [NaN])
    LQOI.set_constraint_primal_solution!(solver, [NaN, NaN])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [NaN])
    LQOI.set_constraint_dual_solution!(solver, [1.0, -1.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    nothing
end