function set_qp1test_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
        termination_status = MOI.Success,
        primal_status = MOI.FeasiblePoint,
        dual_status = MOI.FeasiblePoint,
        variable_primal = [0.571429, 0.428572, 0.857143],
        constraint_primal = [4.0, 1.0],
        variable_dual = [0.0, 0.0, 0.0],
        constraint_dual = [0.714286, 0.857144]
    )
    return nothing
end

function set_qp2test_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
        termination_status = MOI.Success,
        primal_status = MOI.FeasiblePoint,
        dual_status = MOI.FeasiblePoint,
        variable_primal = [0.571429, 0.428572, 0.857143],
        constraint_primal = [4.0, 1.0],
        variable_dual = [0.0, 0.0, 0.0],
        constraint_dual = [0.714286, 0.857144]
    )
    # SOLVE 2
    LQOI.set_solution!(solver,
        termination_status = MOI.Success,
        primal_status = MOI.FeasiblePoint,
        dual_status = MOI.FeasiblePoint,
        variable_primal = [0.571429, 0.428572, 0.857143],
        constraint_primal = [4.0, 1.0],
        variable_dual = [-0.0, -0.0, -0.0],
        constraint_dual = [-1.42857, -1.71429]
    )
    return nothing
end

function set_qp3test_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
        termination_status = MOI.Success,
        primal_status = MOI.FeasiblePoint,
        dual_status = MOI.FeasiblePoint,
        variable_primal = [0.250002, 0.749998],
        constraint_primal = [1.0],
        variable_dual = [6.97141e-6, 4.06505e-8],
        constraint_dual = [2.75]
    )
    # SOLVE 2
    LQOI.set_solution!(solver,
        termination_status = MOI.Success,
        primal_status = MOI.FeasiblePoint,
        dual_status = MOI.FeasiblePoint,
        variable_primal = [1.0, 0.0],
        constraint_primal = [1.0],
        variable_dual = [-0.0, -1.0],
        constraint_dual = [2.0]
    )
    return nothing
end

function set_qcp1test_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
        termination_status = MOI.Success,
        primal_status = MOI.FeasiblePoint,
        dual_status = MOI.FeasiblePoint,
        variable_primal = [0.500003, 1.75],
        constraint_primal = [1.25, 2.25],
        quadratic_primal = [2.0],
        variable_dual = [-0.0, -0.0],
        constraint_dual = [4.48497e-8, 5.48002e-8],
        quadratic_dual = [1.0]
    )
    return nothing
end

function set_qcp2test_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
        termination_status = MOI.Success,
        primal_status = MOI.FeasiblePoint,
        dual_status = MOI.FeasiblePoint,
        variable_primal = [1.41421],
        constraint_primal = Float64[],
        quadratic_primal = [2.0],
        variable_dual = [-0.0],
        constraint_dual = Float64[],
        quadratic_dual = [0.353553]
    )
    return nothing
end

function set_qcp3test_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
        termination_status = MOI.Success,
        primal_status = MOI.FeasiblePoint,
        dual_status = MOI.FeasiblePoint,
        variable_primal = [1.41421],
        constraint_primal = Float64[],
        quadratic_primal = [2.0],
        variable_dual = [0.0],
        constraint_dual = Float64[],
        quadratic_dual = [-0.353553]
    )
    return nothing
end

function set_socp1test_solutions!(solver)
    # SOLVE 1
    LQOI.set_solution!(solver,
        termination_status = MOI.Success,
        primal_status = MOI.FeasiblePoint,
        dual_status = MOI.FeasiblePoint,
        variable_primal = [0.5, 0.5, 0.707107],
        constraint_primal = [1.0],
        quadratic_primal = [-9.55708e-13],
        variable_dual = [-3.72382e-10, -3.72383e-10, 2.4911e-10],
        constraint_dual = [0.707107],
        quadratic_dual = [-1.41421]
    )
    return nothing
end