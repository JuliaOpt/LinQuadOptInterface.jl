function set_linear1test_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
        variable_primal = [1.0, 0.0],
        constraint_primal = [1.0],
        variable_dual = [0.0, 1.0],
        constraint_dual = [-1.0]
    )
    # SOLVE 2
    LQOI.set_solution!(solver,
        variable_primal = [1.0, 0.0],
        constraint_primal = [1.0],
        variable_dual = [0.0, -1.0],
        constraint_dual = [1.0]
    )
    # SOLVE 3
    LQOI.set_solution!(solver,
        variable_primal = [0.0, 0.0, 1.0],
        constraint_primal = [1.0],
        variable_dual = [-1.0, -2.0, 0.0],
        constraint_dual = [2.0]
    )
    # SOLVE 4
    LQOI.set_solution!(solver,
        variable_primal = [-1.0, 0.0, 2.0],
        constraint_primal = [1.0],
        variable_dual = [1.0, 2.0, 0.0],
        constraint_dual = [-2.0]
    )
    # SOLVE 5
    LQOI.set_solution!(solver,
        variable_primal = [1.0, 0.0, 0.0],
        constraint_primal = [1.0],
        variable_dual = [0.0, 1.0, -1.0],
        constraint_dual = [-1.0]
    )
    # SOLVE 6
    LQOI.set_solution!(solver,
        variable_primal = [2.0, 0.0, 0.0],
        constraint_primal = [2.0],
        variable_dual = [0.0, 1.0, -1.0],
        constraint_dual = [-1.0]
    )
    # SOLVE 7
    LQOI.set_solution!(solver,
        variable_primal = [0.0, 2.0, 0.0],
        constraint_primal = [2.0],
        variable_dual = [1.0, 0.0, 2.0],
        constraint_dual = [-2.0]
    )
    # SOLVE 8
    LQOI.set_solution!(solver,
        variable_primal = [1.0, 1.0, 0.0],
        constraint_primal = [2.0, 0.0],
        variable_dual = [0.0, 0.0, -1.5],
        constraint_dual = [1.5, -0.5]
    )
    return
end

function set_linear2test_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
        variable_primal = [1.0, 0.0],
        constraint_primal = [1.0],
        variable_dual = [0.0, 1.0],
        constraint_dual = [-1.0]
    )
    return
end

function set_linear3test_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
        variable_primal = [3.0],
        constraint_primal = [3.0],
        variable_dual = [0.0],
        constraint_dual = [1.0]
    )
    # SOLVE 2
    LQOI.set_solution!(solver,
        variable_primal = [0.0],
        constraint_primal = [0.0],
        variable_dual = [1.0],
        constraint_dual = [0.0]
    )
    return
end

function set_linear4test_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
        variable_primal = [0.0, 0.0],
        constraint_primal = Float64[],
        variable_dual = [1.0, -1.0],
        constraint_dual = Float64[]
    )
    # SOLVE 2
    LQOI.set_solution!(solver,
        variable_primal = [100.0, 0.0],
        constraint_primal = Float64[],
        variable_dual = [1.0, -1.0],
        constraint_dual = Float64[]
    )
    # SOLVE 3
    LQOI.set_solution!(solver,
        variable_primal = [100.0, -100.0],
        constraint_primal = Float64[],
        variable_dual = [1.0, -1.0],
        constraint_dual = Float64[]
    )
    return
end

function set_linear5test_solutions!(solver)

    # SOLVE 1
    LQOI.set_solution!(solver,
        variable_primal = [4/3, 4/3],
        constraint_primal = [4.0, 4.0],
        variable_dual = [-0.0, -0.0],
        constraint_dual = [1/3, 1/3]
    )
    # SOLVE 2
    LQOI.set_solution!(solver,
        variable_primal = [2.0, 0.0],
        constraint_primal = [4.0, 2.0],
        variable_dual = [-0.0, -0.5],
        constraint_dual = [0.5, 0.0]
    )
    # SOLVE 3
    LQOI.set_solution!(solver,
        variable_primal = [4.0, 0.0],
        constraint_primal = [4.0],
        variable_dual = [-0.0, -1.0],
        constraint_dual = [1.0]
    )
    # SOLVE 4
    LQOI.set_solution!(solver,
        variable_primal = [2.0],
        constraint_primal = [4.0],
        variable_dual = [-0.0],
        constraint_dual = [0.5]
    )
    return
end

function set_linear6test_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
        variable_primal = [0.0, 0.0],
        constraint_primal = [0.0, 0.0],
        variable_dual = [0.0, 0.0],
        constraint_dual = [1.0, -1.0]
    )
    # SOLVE 2
    LQOI.set_solution!(solver,
        variable_primal = [100.0, 0.0],
        constraint_primal = [100.0, 0.0],
        variable_dual = [0.0, 0.0],
        constraint_dual = [1.0, -1.0]
    )
    # SOLVE 3
    LQOI.set_solution!(solver,
        variable_primal = [100.0, -100.0],
        constraint_primal = [100.0, -100.0],
        variable_dual = [0.0, 0.0],
        constraint_dual = [1.0, -1.0]
    )
    return
end

function set_linear7test_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
        variable_primal = [-0.0, -0.0],
        constraint_primal = [-0.0, -0.0],
        variable_dual = [0.0, 0.0],
        constraint_dual = [1.0, -1.0]
    )
    # SOLVE 2
    LQOI.set_solution!(solver,
        variable_primal = [100.0, -0.0],
        constraint_primal = [100.0, -0.0],
        variable_dual = [0.0, 0.0],
        constraint_dual = [1.0, -1.0]
    )
    # SOLVE 3
    LQOI.set_solution!(solver,
        variable_primal = [100.0, -100.0],
        constraint_primal = [100.0, -100.0],
        variable_dual = [0.0, 0.0],
        constraint_dual = [1.0, -1.0]
    )
    return
end

function set_linear8atest_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
        termination_status = MOI.Success,
        primal_status = MOI.UnknownResultStatus,
        dual_status = MOI.InfeasibilityCertificate,
        variable_primal = [NaN, NaN],
        constraint_primal = [NaN],
        variable_dual = [NaN, NaN],
        constraint_dual = [-0.5]
    )
    return
end

function set_linear8btest_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
        termination_status = MOI.Success,
        primal_status = MOI.InfeasibilityCertificate,
        dual_status = MOI.UnknownResultStatus,
        variable_primal = [2.0, 1.0],
        constraint_primal = [NaN],
        variable_dual = [NaN, NaN],
        constraint_dual = [NaN]
    )
    return
end

function set_linear8ctest_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
        termination_status = MOI.Success,
        primal_status = MOI.InfeasibilityCertificate,
        dual_status = MOI.UnknownResultStatus,
        variable_primal = [1.0, 1.0],
        constraint_primal = [NaN],
        variable_dual = [NaN, NaN],
        constraint_dual = [NaN]
    )
    return
end

function set_linear9test_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
        variable_primal = [650/11, 400/11],
        constraint_primal = [4.54545, 1000.0, 70000.0],
        variable_dual = [-0.0, -0.0],
        constraint_dual = [0.0, 11.3636, 0.863636]
    )
    return
end

function set_linear10test_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
        variable_primal = [10.0, 0.0],
        constraint_primal = [10.0],
        variable_dual = [-0.0, -0.0],
        constraint_dual = [1.0]
    )
    # SOLVE 2
    LQOI.set_solution!(solver,
        variable_primal = [5.0, 0.0],
        constraint_primal = [5.0],
        variable_dual = [0.0, 0.0],
        constraint_dual = [1.0]
    )
    # SOLVE 3
    LQOI.set_solution!(solver,
        variable_primal = [2.0, 0.0],
        constraint_primal = [2.0],
        variable_dual = [0.0, 0.0],
        constraint_dual = [1.0]
    )
    # SOLVE 4
    LQOI.set_solution!(solver,
        variable_primal = [12.0, 0.0],
        constraint_primal = [12.0],
        variable_dual = [-0.0, -0.0],
        constraint_dual = [1.0]
    )
    return
end

function set_linear11test_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
        variable_primal = [2.0, 0.0],
        constraint_primal = [2.0, 2.0],
        variable_dual = [0.0, 0.0],
        constraint_dual = [0.0, 1.0]
    )
    # SOLVE 2
    LQOI.set_solution!(solver,
        variable_primal = [1.0, 0.0],
        constraint_primal = [1.0, 1.0],
        variable_dual = [0.0, 0.0],
        constraint_dual = [1.0, 0.0]
    )
    return
end

function set_linear12test_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
        termination_status = MOI.Success,
        primal_status = MOI.UnknownResultStatus,
        dual_status = MOI.InfeasibilityCertificate,
        variable_primal = [NaN, NaN],
        constraint_primal = [NaN, NaN],
        variable_dual = [2/3, 0.0],
        constraint_dual = [-1/3, -1.0]
    )
    return
end

function set_linear13test_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
        variable_primal = [0.2, 0.2],
        constraint_primal = [1.0, 0.0],
        variable_dual = [0.0, 0.0],
        constraint_dual = [0.0, 0.0]
    )
    return
end

function set_linear14test_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
        variable_primal = [0.0, 0.5, 1.0],
        constraint_primal = [2.0],
        variable_dual = [-2.0, -0.0, 2.0],
        constraint_dual = [1.0]
    )
    # SOLVE 2
    LQOI.set_solution!(solver,
        variable_primal = [1.0],
        constraint_primal = [2.0],
        variable_dual = [-0.0],
        constraint_dual = [1.0]
    )
    return
end
