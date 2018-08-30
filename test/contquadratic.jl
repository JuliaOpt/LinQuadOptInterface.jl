function set_qp1test_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [0.571429, 0.428572, 0.857143])
    LQOI.set_constraint_primal_solution!(solver, [4.0, 1.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [0.0, 0.0, 0.0])
    LQOI.set_constraint_dual_solution!(solver, [0.714286, 0.857144])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    nothing
end

function set_qp2test_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [0.571429, 0.428572, 0.857143])
    LQOI.set_constraint_primal_solution!(solver, [4.0, 1.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [0.0, 0.0, 0.0])
    LQOI.set_constraint_dual_solution!(solver, [0.714286, 0.857144])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    # SOLVE 2
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [0.571429, 0.428572, 0.857143])
    LQOI.set_constraint_primal_solution!(solver, [4.0, 1.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [-0.0, -0.0, -0.0])
    LQOI.set_constraint_dual_solution!(solver, [-1.42857, -1.71429])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    nothing
end

function set_qp3test_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [0.250002, 0.749998])
    LQOI.set_constraint_primal_solution!(solver, [1.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [6.97141e-6, 4.06505e-8])
    LQOI.set_constraint_dual_solution!(solver, [2.75])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    # SOLVE 2
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [1.0, 0.0])
    LQOI.set_constraint_primal_solution!(solver, [1.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [-0.0, -1.0])
    LQOI.set_constraint_dual_solution!(solver, [2.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    nothing
end

function set_qcp1test_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [0.500003, 1.75])
    LQOI.set_constraint_primal_solution!(solver, [1.25, 2.25])
    LQOI.set_quadratic_primal_solution!(solver, [2.0])

    LQOI.set_variable_dual_solution!(solver, [-0.0, -0.0])
    LQOI.set_constraint_dual_solution!(solver, [4.48497e-8, 5.48002e-8])
    LQOI.set_quadratic_dual_solution!(solver, [1.0])

    nothing
end

function set_qcp2test_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [1.41421])
    LQOI.set_constraint_primal_solution!(solver, Float64[])
    LQOI.set_quadratic_primal_solution!(solver, [2.0])

    LQOI.set_variable_dual_solution!(solver, [-0.0])
    LQOI.set_constraint_dual_solution!(solver, Float64[])
    LQOI.set_quadratic_dual_solution!(solver, [0.353553])

    nothing
end

function set_qcp3test_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [1.41421])
    LQOI.set_constraint_primal_solution!(solver, Float64[])
    LQOI.set_quadratic_primal_solution!(solver, [2.0])

    LQOI.set_variable_dual_solution!(solver, [0.0])
    LQOI.set_constraint_dual_solution!(solver, Float64[])
    LQOI.set_quadratic_dual_solution!(solver, [-0.353553])

    nothing
end

function set_socp1test_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [0.5, 0.5, 0.707107])
    LQOI.set_constraint_primal_solution!(solver, [1.0])
    LQOI.set_quadratic_primal_solution!(solver, [-9.55708e-13])

    LQOI.set_variable_dual_solution!(solver, [-3.72382e-10, -3.72383e-10, 2.4911e-10])
    LQOI.set_constraint_dual_solution!(solver, [0.707107])
    LQOI.set_quadratic_dual_solution!(solver, [-1.41421])

    nothing
end