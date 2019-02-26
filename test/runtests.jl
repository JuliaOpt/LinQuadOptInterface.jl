using Compat.Test, MathOptInterface
using LinQuadOptInterface

const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities
const MOIT = MathOptInterface.Test
const MOIU = MathOptInterface.Utilities
const LQOI = LinQuadOptInterface

# Test that using LinQuadOptInterface passes without error. We have to do this
# here because types can't be created in testsets. If it can be defined without
# error, we're okay. This is most likely to error if fields are added to
# @LinQuadOptimizerBase without proper module prefixes.
using Compat  # For Nothing on v0.6
mutable struct OptimizerTest <: LinQuadOptInterface.LinQuadOptimizer
    LinQuadOptInterface.@LinQuadOptimizerBase
    OptimizerTest() = new()
end

@testset "LinQuadOptInterface" begin
    solver = LQOI.MockLinQuadOptimizer()

    @testset "Printing" begin
        @test sprint(show, solver) == "A LinQuadOptInterface model with " *
            "backend:\nMockLinQuadModel\n    Sense       minimize\n    " *
            "Variables   0\n    Constraints 0\n"
    end

    config = MOIT.TestConfig(solve=false)

    @testset "Unit Tests" begin
        config = MOIT.TestConfig(solve=false)
        MOIT.basic_constraint_tests(solver, config)
        MOIT.unittest(solver, config, [
            "solve_affine_interval",
            "solve_qp_edge_cases",
            "solve_qcp_edge_cases",
            "solve_affine_deletion_edge_cases"
        ])
    end

    @testset "Linear tests" begin
        config = MOIT.TestConfig(solve=false)
        MOIT.contlineartest(solver, config, [
            # partial_start requires VariablePrimalStart to be implemented by the
            # solver.
            "partial_start"
        ])
        config = MOIT.TestConfig()
        include("contlinear.jl")
        set_linear1test_solutions!(solver)
        MOIT.linear1test(solver, config)
        set_linear2test_solutions!(solver)
        MOIT.linear2test(solver, config)
        set_linear3test_solutions!(solver)
        MOIT.linear3test(solver, config)
        set_linear4test_solutions!(solver)
        MOIT.linear4test(solver, config)
        set_linear5test_solutions!(solver)
        MOIT.linear5test(solver, config)
        set_linear6test_solutions!(solver)
        MOIT.linear6test(solver, config)
        set_linear7test_solutions!(solver)
        MOIT.linear7test(solver, config)
        set_linear8atest_solutions!(solver)
        MOIT.linear8atest(solver, config)
        set_linear8btest_solutions!(solver)
        MOIT.linear8btest(solver, config)
        set_linear8ctest_solutions!(solver)
        MOIT.linear8ctest(solver, config)
        set_linear9test_solutions!(solver)
        MOIT.linear9test(solver, config)
        set_linear10test_solutions!(solver)
        MOIT.linear10test(solver, config)
        set_linear11test_solutions!(solver)
        MOIT.linear11test(solver, config)
        set_linear12test_solutions!(solver)
        MOIT.linear12test(solver, MOIT.TestConfig(atol=1e-3,rtol=1e-3))
        set_linear13test_solutions!(solver)
        MOIT.linear13test(solver, config)
        set_linear14test_solutions!(solver)
        MOIT.linear14test(solver, config)
        # set_linear15test_solutions!(solver)
        # MOIT.linear15test(solver, config)
    end

    @testset "Quadratic tests" begin
        config = MOIT.TestConfig(solve=false)
        MOIT.contquadratictest(solver, config)
        include("contquadratic.jl")
        config = MOIT.TestConfig(atol=1e-3, rtol=1e-3, duals=false)
        set_qp1test_solutions!(solver)
        MOIT.qp1test(solver, config)
        set_qp2test_solutions!(solver)
        MOIT.qp2test(solver, config)
        set_qp3test_solutions!(solver)
        MOIT.qp3test(solver, config)
        set_qcp1test_solutions!(solver)
        MOIT.qcp1test(solver, config)
        set_qcp2test_solutions!(solver)
        MOIT.qcp2test(solver, config)
        set_qcp3test_solutions!(solver)
        MOIT.qcp3test(solver, config)
        set_socp1test_solutions!(solver)
        MOIT.socp1test(solver, config)
    end

    @testset "Linear Conic tests" begin
        config = MOIT.TestConfig(solve=false)
        MOIT.lintest(solver, config)
        include("contconic.jl")
        config = MOIT.TestConfig()
        set_lin1test_solutions!(solver)
        MOIT.lin1vtest(solver, config)
        set_lin1test_solutions!(solver)
        MOIT.lin1ftest(solver, config)
        set_lin2test_solutions!(solver)
        MOIT.lin2vtest(solver, config)
        set_lin2test_solutions!(solver)
        MOIT.lin2ftest(solver, config)
        set_lin3test_solutions!(solver)
        MOIT.lin3test(solver, config)
        set_lin3test_solutions!(solver)
        MOIT.lin4test(solver, config)
    end

    @testset "Integer Linear tests" begin
        config = MOIT.TestConfig(solve=false)
        MOIT.intlineartest(solver, config, ["int2"])
        include("intlinear.jl")
        config = MOIT.TestConfig()
        set_knapsacktest_solutions!(solver)
        MOIT.knapsacktest(solver, config)
        set_int1test_solutions!(solver)
        MOIT.int1test(solver, config)
        set_int3test_solutions!(solver)
        MOIT.int3test(solver, config)
    end

    @testset "ModelLike tests" begin
        config = MOIT.TestConfig(solve=false)
        @testset "default_objective_test" begin
            MOIT.default_objective_test(LQOI.MockLinQuadOptimizer())
        end
        @testset "default_status_test" begin
            MOIT.default_status_test(LQOI.MockLinQuadOptimizer())
        end
        @testset "nametest" begin
            MOIT.nametest(LQOI.MockLinQuadOptimizer())
        end
        @testset "validtest" begin
            MOIT.validtest(LQOI.MockLinQuadOptimizer())
        end
        @testset "emptytest" begin
            MOIT.emptytest(LQOI.MockLinQuadOptimizer())
        end
        @testset "orderedindicestest" begin
            MOIT.orderedindicestest(LQOI.MockLinQuadOptimizer())
        end
        @testset "copytest" begin
            MOIT.copytest(LQOI.MockLinQuadOptimizer(),
                          LQOI.MockLinQuadOptimizer())
            MOIT.failcopytestc(LQOI.MockLinQuadOptimizer())
            MOIT.failcopytestia(LQOI.MockLinQuadOptimizer())
            MOIT.failcopytestva(LQOI.MockLinQuadOptimizer())
            MOIT.failcopytestca(LQOI.MockLinQuadOptimizer())
        end
    end
end

MOIU.@model(ModelComplete,
    (MOI.ZeroOne, MOI.Integer),
    (MOI.EqualTo, MOI.GreaterThan, MOI.LessThan, MOI.Interval,
     MOI.Semicontinuous, MOI.Semiinteger),
    (MOI.Reals, MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives,
     MOI.SecondOrderCone, MOI.RotatedSecondOrderCone, MOI.GeometricMeanCone,
     MOI.ExponentialCone, MOI.DualExponentialCone,
     MOI.PositiveSemidefiniteConeTriangle, MOI.PositiveSemidefiniteConeSquare,
     MOI.RootDetConeTriangle, MOI.RootDetConeSquare, MOI.LogDetConeTriangle,
     MOI.LogDetConeSquare),
    (MOI.PowerCone, MOI.DualPowerCone, MOI.SOS1, MOI.SOS2),
    (MOI.SingleVariable,),
    (MOI.ScalarAffineFunction, MOI.ScalarQuadraticFunction),
    (MOI.VectorOfVariables,),
    (MOI.VectorAffineFunction, MOI.VectorQuadraticFunction))
@testset "Copy from/to @Model" begin
    MOIT.copytest(LQOI.MockLinQuadOptimizer(), ModelComplete{Float64}())
    MOIT.copytest(ModelComplete{Float64}(), LQOI.MockLinQuadOptimizer())
end

@testset "Issue #52" begin
    @testset "scalaraffine" begin
        model = LQOI.MockLinQuadOptimizer()
        x = MOI.add_variable(model)
        f1 =  MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0], [x]), 1.5)
        s1 = MOI.LessThan(2.0)
        c = MOI.add_constraint(model, f1, s1)
        f2 = MOI.get(model, MOI.ConstraintFunction(), c)
        @test f1 ≈ f2
        s2 = MOI.get(model, MOI.ConstraintSet(), c)
        @test typeof(s1) == typeof(s2)
        @test s1.upper == s2.upper
    end
    @testset "vectoraffine" begin
        model = LQOI.MockLinQuadOptimizer()
        x = MOI.add_variable(model)
        f1 =  MOI.VectorAffineFunction(
            MOI.VectorAffineTerm.([1], MOI.ScalarAffineTerm.([1.0], [x])), [1.5])
        s1 = MOI.Zeros(1)
        c = MOI.add_constraint(model, f1, s1)
        f2 = MOI.get(model, MOI.ConstraintFunction(), c)
        @test f1 ≈ f2
        s2 = MOI.get(model, MOI.ConstraintSet(), c)
        @test typeof(s1) == typeof(s2)
        @test s1.dimension == s2.dimension
    end
end

@testset "Issue #54" begin
    @testset "scalar, one-sided" begin
        model = LQOI.MockLinQuadOptimizer()
        x = MOI.add_variable(model)
        f =  MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0], [x]), 1.5)
        s = MOI.LessThan(2.0)
        c = MOI.add_constraint(model, f, s)
        # Change the constraint set, and verify that we get the same set
        # when we retrieve it:
        s2 = MOI.LessThan(1.0)
        MOI.set(model, MOI.ConstraintSet(), c, s2)
        s3 = MOI.get(model, MOI.ConstraintSet(), c)
        @test typeof(s2) == typeof(s3)
        @test s2.upper == s3.upper
    end
    @testset "scalar, interval" begin
        model = LQOI.MockLinQuadOptimizer()
        x = MOI.add_variable(model)
        f =  MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0], [x]), 1.5)
        s = MOI.Interval(2.0, 3.0)
        c = MOI.add_constraint(model, f, s)
        # Change the constraint set, and verify that we get the same set
        # when we retrieve it:
        s2 = MOI.Interval(1.0, 2.0)
        MOI.set(model, MOI.ConstraintSet(), c, s2)
        s3 = MOI.get(model, MOI.ConstraintSet(), c)
        @test typeof(s2) == typeof(s3)
        @test s2.lower == s3.lower
        @test s2.upper == s3.upper
    end
end

@testset "Issue #32" begin
    check_set(actual::MOI.GreaterThan, expected::MOI.GreaterThan) =
        @test actual.lower == expected.lower
    check_set(actual::MOI.LessThan, expected::MOI.LessThan) =
        @test actual.upper == expected.upper
    function check_set(actual::MOI.Interval, expected::MOI.Interval)
        @test actual.lower == expected.lower
        @test actual.lower == expected.lower
    end
    function check_set(actual::S, expected::S) where {S <: Union{MOI.Nonnegatives, MOI.Nonpositives, MOI.Zeros}}
        @test actual.dimension == expected.dimension
    end

    function check_row(model, constraint_index, expected_function, expected_set)
        actual_function = MOI.get(model, MOI.ConstraintFunction(), constraint_index)
        @test actual_function ≈ expected_function
        actual_set = MOI.get(model, MOI.ConstraintSet(), constraint_index)
        check_set(actual_set, expected_set)
    end

    @testset "Scalar affine function" begin
        @testset "Single variable" begin
            model = LQOI.MockLinQuadOptimizer()
            x = MOI.add_variable(model)

            # Add the constraint and verify that we can retrieve it
            f =  MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([2.0], [x]), 0.1)
            s = MOI.GreaterThan(0.5)
            c = MOI.add_constraint(model, f, s)
            check_row(model, c, f, s)

            # Replace the constraint function with itself and verify that
            # the problem is unchanged
            MOI.set(model, MOI.ConstraintFunction(), c, f)
            check_row(model, c, f, s)

            # Replace the constraint function with a new function, verify
            # that the replacement occurred and that the set is unchanged
            f2 = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0], [x]), 0.5)
            MOI.set(model, MOI.ConstraintFunction(), c, f2)
            check_row(model, c, f2, s)

            # Replace the constraint function with a new function which is not
            # in canonical form:
            f3 = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0, 0.5], [x]), 0.5)
            MOI.set(model, MOI.ConstraintFunction(), c, f3)
            check_row(model, c, MOIU.canonical(f3), s)

            # Replace the constraint function with a new function whose sparsity
            # pattern is different:
            f4 = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Float64}[], 2.0)
            MOI.set(model, MOI.ConstraintFunction(), c, f4)
            check_row(model, c, f4, s)
        end

        @testset "Multiple variables" begin
            model = LQOI.MockLinQuadOptimizer()
            x = MOI.add_variable(model)
            y = MOI.add_variable(model)

            # Add the constraint and verify that we can retrieve it
            f =  MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0, 2.0], [x, y]), 0.1)
            s = MOI.GreaterThan(0.5)
            c = MOI.add_constraint(model, f, s)
            check_row(model, c, f, s)

            # Replace the constraint function with itself and verify that
            # the problem is unchanged
            MOI.set(model, MOI.ConstraintFunction(), c, f)
            check_row(model, c, f, s)

            # Replace the constraint function with a new function, verify
            # that the replacement occurred and that the set is unchanged
            f2 = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([0.5, 1.0], [x, y]), 0.5)
            MOI.set(model, MOI.ConstraintFunction(), c, f2)
            check_row(model, c, f2, s)

            # Replace the constraint function with a new function which is not
            # in canonical form:
            f3 = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([-0.1, 0.2], [y, x]), 0.5)
            MOI.set(model, MOI.ConstraintFunction(), c, f3)
            check_row(model, c, MOIU.canonical(f3), s)

            # Replace the constraint function with a new function whose sparsity
            # pattern is different:
            f4 = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([-2.0], [y]), 2.0)
            MOI.set(model, MOI.ConstraintFunction(), c, f4)
            check_row(model, c, f4, s)
        end
    end

    @testset "Vector affine function" begin
        model = LQOI.MockLinQuadOptimizer()
        x = MOI.add_variable(model)
        y = MOI.add_variable(model)

        # Add the constraint and verify that we can retrieve it
        f = MOI.VectorAffineFunction(
            [MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(1.0, x)),
             MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(2.0, y)),
             MOI.VectorAffineTerm(2, MOI.ScalarAffineTerm(-1.0, x))],
            [0.1, 0.2])
        s = MOI.Nonpositives(2)
        c = MOI.add_constraint(model, f, s)
        check_row(model, c, f, s)

        # Replace the constraint function with itself and verify that
        # the problem is unchanged
        MOI.set(model, MOI.ConstraintFunction(), c, f)
        check_row(model, c, f, s)


        # Replace the constraint function with a new function, verify
        # that the replacement occurred and that the set is unchanged
        f2 = MOI.VectorAffineFunction(
            [MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(3.0, x)),
             MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(-2.0, y)),
             MOI.VectorAffineTerm(2, MOI.ScalarAffineTerm(-1.5, x))],
            [0.5, 2.0])
        MOI.set(model, MOI.ConstraintFunction(), c, f2)
        check_row(model, c, f2, s)

        # Replace the constraint function with a new function which is not
        # in canonical form:
        f3 = MOI.VectorAffineFunction(
            [MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(2.0, x)),
             MOI.VectorAffineTerm(2, MOI.ScalarAffineTerm(-1.5, x)),
             MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(-2.5, y))],
            [0.1, 0.2])
        MOI.set(model, MOI.ConstraintFunction(), c, f3)
        check_row(model, c, MOIU.canonical(f3), s)

        # Replace the constraint function with a new function whose sparsity
        # pattern is different:
        f4 = MOI.VectorAffineFunction(
            [MOI.VectorAffineTerm(2, MOI.ScalarAffineTerm(-1.5, y)),
             MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(2.0, x))],
            [0.1, 0.2])
        MOI.set(model, MOI.ConstraintFunction(), c, f4)
        check_row(model, c, f4, s)
    end
end

@testset "Issue #67" begin
    model = LQOI.MockLinQuadOptimizer()
    function objective_type_test(f::MOI.AbstractScalarFunction)
        MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
        @test MOI.get(model, MOI.ObjectiveFunctionType()) == typeof(f)
    end
    x = MOI.add_variable(model)
    f = MOI.SingleVariable(x)
    objective_type_test(f)
    objective_type_test(convert(MOI.ScalarAffineFunction{Float64}, f))
    objective_type_test(convert(MOI.ScalarQuadraticFunction{Float64}, f))
end

@testset "Conflicting SingleVariable constraints" begin
    @testset "ZeroOne" begin
        model = LQOI.MockLinQuadOptimizer()
        x = MOI.add_variable(model)
        MOI.add_constraint(model, MOI.SingleVariable(x), MOI.ZeroOne())
        @test_throws Exception MOI.add_constraint(model, MOI.SingleVariable(x), MOI.Interval(2.0, 3.0))
        @test_throws Exception MOI.add_constraint(model, MOI.SingleVariable(x), MOI.Integer())
        @test_throws Exception MOI.add_constraint(model, MOI.SingleVariable(x), MOI.Semiinteger(2.0, 3.0))
        @test_throws Exception MOI.add_constraint(model, MOI.SingleVariable(x), MOI.Semicontinuous(2.0, 3.0))
    end
    @testset "Integer" begin
        model = LQOI.MockLinQuadOptimizer()
        x = MOI.add_variable(model)
        MOI.add_constraint(model, MOI.SingleVariable(x), MOI.Integer())
        MOI.add_constraint(model, MOI.SingleVariable(x), MOI.Interval(2.0, 3.0))
        @test_throws Exception MOI.add_constraint(model, MOI.SingleVariable(x), MOI.ZeroOne())
        @test_throws Exception MOI.add_constraint(model, MOI.SingleVariable(x), MOI.Semiinteger(2.0, 3.0))
        @test_throws Exception MOI.add_constraint(model, MOI.SingleVariable(x), MOI.Semicontinuous(2.0, 3.0))
    end
    @testset "Semiinteger" begin
        model = LQOI.MockLinQuadOptimizer()
        x = MOI.add_variable(model)
        MOI.add_constraint(model, MOI.SingleVariable(x), MOI.Semiinteger(2.0, 3.0))
        @test_throws Exception MOI.add_constraint(model, MOI.SingleVariable(x), MOI.Interval(2.0, 3.0))
        @test_throws Exception MOI.add_constraint(model, MOI.SingleVariable(x), MOI.ZeroOne())
        @test_throws Exception MOI.add_constraint(model, MOI.SingleVariable(x), MOI.Integer())
        @test_throws Exception MOI.add_constraint(model, MOI.SingleVariable(x), MOI.Semicontinuous(2.0, 3.0))
    end
    @testset "Semicontinuous" begin
        model = LQOI.MockLinQuadOptimizer()
        x = MOI.add_variable(model)
        MOI.add_constraint(model, MOI.SingleVariable(x), MOI.Semicontinuous(2.0, 3.0))
        @test_throws Exception MOI.add_constraint(model, MOI.SingleVariable(x), MOI.Interval(2.0, 3.0))
        @test_throws Exception MOI.add_constraint(model, MOI.SingleVariable(x), MOI.ZeroOne())
        @test_throws Exception MOI.add_constraint(model, MOI.SingleVariable(x), MOI.Integer())
        @test_throws Exception MOI.add_constraint(model, MOI.SingleVariable(x), MOI.Semiinteger(2.0, 3.0))
    end
end

@testset "SemiXXX variables" begin
    @testset "Semiinteger" begin
        model = LQOI.MockLinQuadOptimizer()
        x = MOI.add_variable(model)
        c = MOI.add_constraint(
            model, MOI.SingleVariable(x), MOI.Semiinteger(1.0, 4.5))
        @test model.inner.vartype[1] == Cchar('N')
        @test MOI.get(model, MOI.ConstraintSet(), c) == MOI.Semiinteger(1.0, 4.5)
        @test MOI.get(model, MOI.ConstraintFunction(), c) == MOI.SingleVariable(x)
        MOI.delete(model, c)
        @test model.inner.vartype[1] == Cchar('C')
    end
    @testset "Semicontinuous" begin
        model = LQOI.MockLinQuadOptimizer()
        x = MOI.add_variable(model)
        c = MOI.add_constraint(
            model, MOI.SingleVariable(x), MOI.Semicontinuous(1.0, 4.0))
        @test model.inner.vartype[1] == Cchar('S')
        @test MOI.get(model, MOI.ConstraintSet(), c) == MOI.Semicontinuous(1.0, 4.0)
        @test MOI.get(model, MOI.ConstraintFunction(), c) == MOI.SingleVariable(x)
        MOI.delete(model, c)
        @test model.inner.vartype[1] == Cchar('C')
    end
end

@testset "VariablePrimalStart" begin
    model = LQOI.MockLinQuadOptimizer()
    x = MOI.add_variable(model)
    @test_throws MOI.UnsupportedAttribute MOI.set(model, MOI.VariablePrimalStart(), x, 1.0)
end
