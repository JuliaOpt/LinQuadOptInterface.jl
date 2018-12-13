function set_knapsacktest_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
        termination_status = MOI.Optimal,
        primal_status = MOI.FeasiblePoint,
        dual_status = MOI.NoSolution,
        variable_primal = [1.0, 0.0, 0.0, 1.0, 1.0],
        constraint_primal = [9.0]
    )
    return
end

function set_int1test_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
        termination_status = MOI.Optimal,
        primal_status = MOI.FeasiblePoint,
        dual_status = MOI.NoSolution,
        variable_primal = [4.0, 5.0, 1.0],
        constraint_primal = [10.0, 15.0]
    )
    return
end

function set_int3test_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
        termination_status = MOI.Optimal,
        primal_status = MOI.FeasiblePoint,
        dual_status = MOI.NoSolution,
        variable_primal = [1.0, -0.0, -0.0, -0.0, -0.0, 1.0, -0.0, -0.0, -0.0, -0.0, -0.0],
        constraint_primal = [0.9875]
    )
    # SOLVE 2
    LQOI.set_solution!(solver,
        termination_status = MOI.Optimal,
        primal_status = MOI.FeasiblePoint,
        dual_status = MOI.NoSolution,
        variable_primal = [1.0, -0.0, -0.0, -0.0, -0.0, 1.0, -0.0, -0.0, -0.0, -0.0, -0.0],
        constraint_primal = [0.9875]
    )
    return
end
