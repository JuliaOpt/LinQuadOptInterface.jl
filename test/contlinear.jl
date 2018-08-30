function set_linear1test_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [1.0, 0.0])
    LQOI.set_constraint_primal_solution!(solver, [1.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [0.0, 1.0])
    LQOI.set_constraint_dual_solution!(solver, [-1.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    # SOLVE 2
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [1.0, 0.0])
    LQOI.set_constraint_primal_solution!(solver, [1.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [0.0, -1.0])
    LQOI.set_constraint_dual_solution!(solver, [1.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    # SOLVE 3
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [0.0, 0.0, 1.0])
    LQOI.set_constraint_primal_solution!(solver, [1.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [-1.0, -2.0, 0.0])
    LQOI.set_constraint_dual_solution!(solver, [2.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    # SOLVE 4
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)
    
    LQOI.set_variable_primal_solution!(solver, [-1.0, 0.0, 2.0])
    LQOI.set_constraint_primal_solution!(solver, [1.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [1.0, 2.0, 0.0])
    LQOI.set_constraint_dual_solution!(solver, [-2.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    # SOLVE 5
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [1.0, 0.0, 0.0])
    LQOI.set_constraint_primal_solution!(solver, [1.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [0.0, 1.0, -1.0])
    LQOI.set_constraint_dual_solution!(solver, [-1.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    # SOLVE 6
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [2.0, 0.0, 0.0])
    LQOI.set_constraint_primal_solution!(solver, [2.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [0.0, 1.0, -1.0])
    LQOI.set_constraint_dual_solution!(solver, [-1.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    # SOLVE 7
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [0.0, 2.0, 0.0])
    LQOI.set_constraint_primal_solution!(solver, [2.0])
    LQOI.set_quadratic_primal_solution!(solver,  Float64[])

    LQOI.set_variable_dual_solution!(solver, [1.0, 0.0, 2.0])
    LQOI.set_constraint_dual_solution!(solver, [-2.0])
    LQOI.set_quadratic_dual_solution!(solver,  Float64[])

    # SOLVE 8
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [1.0, 1.0, 0.0])
    LQOI.set_constraint_primal_solution!(solver, [2.0, 0.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [0.0, 0.0, -1.5])
    LQOI.set_constraint_dual_solution!(solver, [1.5, -0.5])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    nothing
end

function set_linear2test_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [1.0, 0.0])
    LQOI.set_constraint_primal_solution!(solver, [1.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [0.0, 1.0])
    LQOI.set_constraint_dual_solution!(solver, [-1.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    nothing
end

function set_linear3test_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [3.0])
    LQOI.set_constraint_primal_solution!(solver, [3.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [0.0])
    LQOI.set_constraint_dual_solution!(solver, [1.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    # SOLVE 2
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [0.0])
    LQOI.set_constraint_primal_solution!(solver, [0.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [1.0])
    LQOI.set_constraint_dual_solution!(solver, [0.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    nothing
end

function set_linear4test_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [0.0, 0.0])
    LQOI.set_constraint_primal_solution!(solver, Float64[])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [1.0, -1.0])
    LQOI.set_constraint_dual_solution!(solver, Float64[])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    # SOLVE 2
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [100.0, 0.0])
    LQOI.set_constraint_primal_solution!(solver, Float64[])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [1.0, -1.0])
    LQOI.set_constraint_dual_solution!(solver, Float64[])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    # SOLVE 3
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [100.0, -100.0])
    LQOI.set_constraint_primal_solution!(solver, Float64[])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [1.0, -1.0])
    LQOI.set_constraint_dual_solution!(solver, Float64[])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    nothing
end

function set_linear5test_solutions!(solver)

    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [4/3, 4/3])
    LQOI.set_constraint_primal_solution!(solver, [4.0, 4.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [-0.0, -0.0])
    LQOI.set_constraint_dual_solution!(solver, [1/3, 1/3])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    # SOLVE 2
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [2.0, 0.0])
    LQOI.set_constraint_primal_solution!(solver, [4.0, 2.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [-0.0, -0.5])
    LQOI.set_constraint_dual_solution!(solver, [0.5, 0.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    # SOLVE 3
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [4.0, 0.0])
    LQOI.set_constraint_primal_solution!(solver, [4.0])
    LQOI.set_quadratic_primal_solution!(solver,  Float64[])

    LQOI.set_variable_dual_solution!(solver, [-0.0, -1.0])
    LQOI.set_constraint_dual_solution!(solver, [1.0])
    LQOI.set_quadratic_dual_solution!(solver,  Float64[])

    # SOLVE 4
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [2.0])
    LQOI.set_constraint_primal_solution!(solver, [4.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [-0.0])
    LQOI.set_constraint_dual_solution!(solver, [0.5])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    nothing
end

function set_linear6test_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [0.0, 0.0])
    LQOI.set_constraint_primal_solution!(solver, [0.0, 0.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [0.0, 0.0])
    LQOI.set_constraint_dual_solution!(solver, [1.0, -1.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    # SOLVE 2
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [100.0, 0.0])
    LQOI.set_constraint_primal_solution!(solver, [100.0, 0.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [0.0, 0.0])
    LQOI.set_constraint_dual_solution!(solver, [1.0, -1.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    # SOLVE 3
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [100.0, -100.0])
    LQOI.set_constraint_primal_solution!(solver, [100.0, -100.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [0.0, 0.0])
    LQOI.set_constraint_dual_solution!(solver, [1.0, -1.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    nothing
end

function set_linear7test_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [-0.0, -0.0])
    LQOI.set_constraint_primal_solution!(solver, [-0.0, -0.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [0.0, 0.0])
    LQOI.set_constraint_dual_solution!(solver, [1.0, -1.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    # SOLVE 2
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [100.0, -0.0])
    LQOI.set_constraint_primal_solution!(solver, [100.0, -0.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [0.0, 0.0])
    LQOI.set_constraint_dual_solution!(solver, [1.0, -1.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    # SOLVE 3
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [100.0, -100.0])
    LQOI.set_constraint_primal_solution!(solver, [100.0, -100.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [0.0, 0.0])
    LQOI.set_constraint_dual_solution!(solver, [1.0, -1.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])
    nothing
end

function set_linear8atest_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.UnknownResultStatus)
    LQOI.set_dual_status!(solver, MOI.InfeasibilityCertificate)

    LQOI.set_variable_primal_solution!(solver, [NaN, NaN])
    LQOI.set_constraint_primal_solution!(solver, [NaN])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [NaN, NaN])
    LQOI.set_constraint_dual_solution!(solver, [-0.5])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    nothing
end

function set_linear8btest_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.InfeasibilityCertificate)
    LQOI.set_dual_status!(solver, MOI.UnknownResultStatus)

    LQOI.set_variable_primal_solution!(solver, [2.0, 1.0])
    LQOI.set_constraint_primal_solution!(solver, [NaN])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [NaN, NaN])
    LQOI.set_constraint_dual_solution!(solver, [NaN])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    nothing
end

function set_linear8ctest_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.InfeasibilityCertificate)
    LQOI.set_dual_status!(solver, MOI.UnknownResultStatus)

    LQOI.set_variable_primal_solution!(solver, [1.0, 1.0])
    LQOI.set_constraint_primal_solution!(solver, [NaN])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [NaN, NaN])
    LQOI.set_constraint_dual_solution!(solver, [NaN])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    nothing
end

function set_linear9test_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [650/11, 400/11])
    LQOI.set_constraint_primal_solution!(solver, [4.54545, 1000.0, 70000.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [-0.0, -0.0])
    LQOI.set_constraint_dual_solution!(solver, [0.0, 11.3636, 0.863636])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    nothing
end

function set_linear10test_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [10.0, 0.0])
    LQOI.set_constraint_primal_solution!(solver, [10.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [-0.0, -0.0])
    LQOI.set_constraint_dual_solution!(solver, [1.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    # SOLVE 2
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [5.0, 0.0])
    LQOI.set_constraint_primal_solution!(solver, [5.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [0.0, 0.0])
    LQOI.set_constraint_dual_solution!(solver, [1.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    # SOLVE 3
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [2.0, 0.0])
    LQOI.set_constraint_primal_solution!(solver, [2.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [0.0, 0.0])
    LQOI.set_constraint_dual_solution!(solver, [1.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    # SOLVE 4
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [12.0, 0.0])
    LQOI.set_constraint_primal_solution!(solver, [12.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [-0.0, -0.0])
    LQOI.set_constraint_dual_solution!(solver, [1.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])
    nothing
end

function set_linear11test_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [2.0, 0.0])
    LQOI.set_constraint_primal_solution!(solver, [2.0, 2.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [0.0, 0.0])
    LQOI.set_constraint_dual_solution!(solver, [0.0, 1.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    # SOLVE 2
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [1.0, 0.0])
    LQOI.set_constraint_primal_solution!(solver, [1.0, 1.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [0.0, 0.0])
    LQOI.set_constraint_dual_solution!(solver, [1.0, 0.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    nothing
end

function set_linear12test_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.UnknownResultStatus)
    LQOI.set_dual_status!(solver, MOI.InfeasibilityCertificate)

    LQOI.set_variable_primal_solution!(solver, [NaN, NaN])
    LQOI.set_constraint_primal_solution!(solver, [NaN, NaN])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [2/3, 0.0])
    LQOI.set_constraint_dual_solution!(solver, [-1/3, -1.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    nothing
end

function set_linear13test_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [0.2, 0.2])
    LQOI.set_constraint_primal_solution!(solver, [1.0, 0.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [0.0, 0.0])
    LQOI.set_constraint_dual_solution!(solver, [0.0, 0.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    nothing
end

function set_linear14test_solutions!(solver)
    # SOLVE 1
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [0.0, 0.5, 1.0])
    LQOI.set_constraint_primal_solution!(solver, [2.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [-2.0, -0.0, 2.0])
    LQOI.set_constraint_dual_solution!(solver, [1.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])

    # SOLVE 2
    LQOI.set_termination_status!(solver, MOI.Success)
    LQOI.set_primal_status!(solver, MOI.FeasiblePoint)
    LQOI.set_dual_status!(solver, MOI.FeasiblePoint)

    LQOI.set_variable_primal_solution!(solver, [1.0])
    LQOI.set_constraint_primal_solution!(solver, [2.0])
    LQOI.set_quadratic_primal_solution!(solver, Float64[])

    LQOI.set_variable_dual_solution!(solver, [-0.0])
    LQOI.set_constraint_dual_solution!(solver, [1.0])
    LQOI.set_quadratic_dual_solution!(solver, Float64[])
    nothing
end
