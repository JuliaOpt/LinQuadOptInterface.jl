using Compat.Test, MathOptInterface
using LinQuadOptInterface

const MOI = MathOptInterface
const MOIT = MathOptInterface.Test
const LQOI = LinQuadOptInterface

@testset "LinQuadOptInterface" begin
    solver = LQOI.MockLinQuadOptimizer()
    # TODO(@joaquim): test with solve=true
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
        MOIT.contlineartest(solver, config)
        config = MOIT.TestConfig(solve=true)
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
        # set_linear9test_solutions!(solver)
        # MOIT.linear9test(solver, config)
        # set_linear10test_solutions!(solver)
        # MOIT.linear10test(solver, config)
        # set_linear11test_solutions!(solver)
        # MOIT.linear11test(solver, config)
        # set_linear12test_solutions!(solver)
        # MOIT.linear12test(solver, config)
        # set_linear13test_solutions!(solver)
        # MOIT.linear13test(solver, config)
        # set_linear14test_solutions!(solver)
        # MOIT.linear14test(solver, config)
        # set_linear15test_solutions!(solver)
        # MOIT.linear15test(solver, config)
    end

    @testset "Quadratic tests" begin
        config = MOIT.TestConfig(solve=false)
        MOIT.contquadratictest(solver, config)
    end

    @testset "Linear Conic tests" begin
        config = MOIT.TestConfig(solve=false)
        MOIT.lintest(solver, config)
    end

    @testset "Integer Linear tests" begin
        config = MOIT.TestConfig(solve=false)
        MOIT.intlineartest(solver, config, ["int2"])
    end

    @testset "ModelLike tests" begin
        config = MOIT.TestConfig(solve=false)
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
        end
    end
end

@testset "Issue #52" begin
    @testset "scalaraffine" begin
        model = LQOI.MockLinQuadOptimizer()
        x = MOI.addvariable!(model)
        f1 =  MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0], [x]), 1.5)
        s1 = MOI.LessThan(2.0)
        c = MOI.addconstraint!(model, f1, s1)
        f2 = MOI.get(model, MOI.ConstraintFunction(), c)
        @test f1 ≈ f2
        s2 = MOI.get(model, MOI.ConstraintSet(), c)
        @test typeof(s1) == typeof(s2)
        @test s1.upper == s2.upper
    end
    @testset "vectoraffine" begin
        model = LQOI.MockLinQuadOptimizer()
        x = MOI.addvariable!(model)
        f1 =  MOI.VectorAffineFunction(
            MOI.VectorAffineTerm.([1], MOI.ScalarAffineTerm.([1.0], [x])), [1.5])
        s1 = MOI.Zeros(1)
        c = MOI.addconstraint!(model, f1, s1)
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
        x = MOI.addvariable!(model)
        f =  MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0], [x]), 1.5)
        s = MOI.LessThan(2.0)
        c = MOI.addconstraint!(model, f, s)
        # Change the constraint set, and verify that we get the same set
        # when we retrieve it:
        s2 = MOI.LessThan(1.0)
        MOI.set!(model, MOI.ConstraintSet(), c, s2)
        s3 = MOI.get(model, MOI.ConstraintSet(), c)
        @test typeof(s2) == typeof(s3)
        @test s2.upper == s3.upper
    end
    @testset "scalar, interval" begin
        model = LQOI.MockLinQuadOptimizer()
        x = MOI.addvariable!(model)
        f =  MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0], [x]), 1.5)
        s = MOI.Interval(2.0, 3.0)
        c = MOI.addconstraint!(model, f, s)
        # Change the constraint set, and verify that we get the same set
        # when we retrieve it:
        s2 = MOI.Interval(1.0, 2.0)
        MOI.set!(model, MOI.ConstraintSet(), c, s2)
        s3 = MOI.get(model, MOI.ConstraintSet(), c)
        @test typeof(s2) == typeof(s3)
        @test s2.lower == s3.lower
        @test s2.upper == s3.upper
    end
end
