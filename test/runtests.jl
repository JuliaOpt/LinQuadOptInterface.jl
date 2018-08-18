using Compat.Test, MathOptInterface
using LinQuadOptInterface

const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities
const MOIT = MathOptInterface.Test
const LQOI = LinQuadOptInterface

@testset "LinQuadOptInterface" begin
    solver = LQOI.MockLinQuadOptimizer()
    # TODO(@joaquim): test with solve=true
    config = MOIT.TestConfig(solve=false)

    @testset "Unit Tests" begin
        MOIT.basic_constraint_tests(solver, config)
        MOIT.unittest(solver, config, [
            "solve_affine_interval",
            "solve_qp_edge_cases",
            "solve_qcp_edge_cases",
            "solve_affine_deletion_edge_cases"
        ])
    end

    @testset "Linear tests" begin
        MOIT.contlineartest(solver, config)
    end

    @testset "Quadratic tests" begin
        MOIT.contquadratictest(solver, config)
    end

    @testset "Linear Conic tests" begin
        MOIT.lintest(solver, config)
    end

    @testset "Integer Linear tests" begin
        MOIT.intlineartest(solver, config, ["int2"])
    end

    @testset "ModelLike tests" begin
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

@testset "Issue #32" begin
    @testset "Simple scalar affine function replacement" begin
        model = LQOI.MockLinQuadOptimizer()
        x = MOI.addvariable!(model)

        # Add the constraint and verify that we can retrieve it
        f =  MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([2.0], [x]), 0.1)
        s = MOI.GreaterThan(0.5)
        c = MOI.addconstraint!(model, f, s)
        @test MOI.get(model, MOI.ConstraintFunction(), c) ≈ f
        s2 = MOI.get(model, MOI.ConstraintSet(), c)
        @test typeof(s2) == typeof(s)
        @test s2.lower == s.lower

        # Replace the constraint function with itself and verify that
        # the problem is unchanged
        MOI.set!(model, MOI.ConstraintFunction(), c, f)
        @test MOI.get(model, MOI.ConstraintFunction(), c) ≈ f
        s2 = MOI.get(model, MOI.ConstraintSet(), c)
        @test typeof(s2) == typeof(s)
        @test s2.lower == s.lower

        # Replace the constraint function with a new function, verify
        # that the replacement occurred and that the set is unchanged
        f2 = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0], [x]), 0.5)
        MOI.set!(model, MOI.ConstraintFunction(), c, f2)
        @test MOI.get(model, MOI.ConstraintFunction(), c) ≈ f2
        s2 = MOI.get(model, MOI.ConstraintSet(), c)
        @test typeof(s2) == typeof(s)
        @test s2.lower == s.lower

        # Replace the constraint function with a new function which is not
        # in canonical form:
        f3 = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([1.0, 0.5], [x]), 0.5)
        MOI.set!(model, MOI.ConstraintFunction(), c, f3)
        @test MOI.get(model, MOI.ConstraintFunction(), c) ≈ MOIU.canonical(f3)
        s2 = MOI.get(model, MOI.ConstraintSet(), c)
        @test typeof(s2) == typeof(s)
        @test s2.lower == s.lower


        # Replace the constraint function with a new function whose sparsity
        # pattern is different:
        f4 = MOI.ScalarAffineFunction(MOI.ScalarAffineTerm{Float64}[], 2.0)
        MOI.set!(model, MOI.ConstraintFunction(), c, f4)
        @test MOI.get(model, MOI.ConstraintFunction(), c) ≈ f4
        s2 = MOI.get(model, MOI.ConstraintSet(), c)
        @test typeof(s2) == typeof(s)
        @test s2.lower == s.lower
    end
end




