function set_lin1test_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
        variable_primal = [1.0, 0.0, 2.0],
        constraint_primal = [1.0, 0.0, 2.0, 3.0, 2.0],
        variable_dual = [0.0, 0.0, 0.0],
        constraint_dual = [-0.0, 2.0, -0.0, -3.0, -1.0]
    )
    return
end

function set_lin2test_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
        variable_primal = [-4.0, -3.0, 16.0, 0.0],
        constraint_primal = [-4.0, -3.0, 12.0, -3.0, 16.0, 0.0],
        variable_dual = [0.0, 0.0, 0.0, 0.0],
        constraint_dual = [7.0, 2.0, -4.0, -0.0, -0.0, 7.0]
    )
    return
end

function set_lin3test_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
        termination_status = MOI.Infeasible,
        primal_status = MOI.NoSolution,
        dual_status = MOI.InfeasibilityCertificate,
        variable_primal = [NaN],
        constraint_primal = [NaN, NaN],
        variable_dual = [NaN],
        constraint_dual = [1.0, -1.0]
    )
    return
end
