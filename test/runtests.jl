using Compat.Test, MathOptInterface
using LinQuadOptInterface

const MOI = MathOptInterface
const MOIT = MathOptInterface.Test
const LQOI = LinQuadOptInterface

@testset "LinQuadOptInterface" begin
    @testset "Unit Tests" begin
        config = MOIT.TestConfig(solve = false)
        solver = LQOI.MockLinQuadOptimizer()
        MOIT.basic_constraint_tests(solver, config;
            exclude = [
                (MOI.SingleVariable, MOI.EqualTo{Float64}),
                (MOI.SingleVariable, MOI.Integer),
                (MOI.SingleVariable, MOI.LessThan{Float64}),
                (MOI.SingleVariable, MOI.Interval{Float64}),
                (MOI.SingleVariable, MOI.GreaterThan{Float64})
            ]
        )
        MOIT.unittest(solver, config, [
            "solve_affine_interval",
            "solve_qp_edge_cases",
            "solve_qcp_edge_cases",
            "solve_affine_deletion_edge_cases"
        ])
    end

    @testset "Linear tests" begin
        linconfig = MOIT.TestConfig(solve = false)
        solver = LQOI.MockLinQuadOptimizer()
        MOIT.contlineartest(solver , linconfig)
    end

    @testset "Quadratic tests" begin
        quadconfig = MOIT.TestConfig(solve=false)
        solver = LQOI.MockLinQuadOptimizer()
        MOIT.contquadratictest(solver, quadconfig)
    end

    @testset "Linear Conic tests" begin
        linconfig = MOIT.TestConfig(solve=false)
        solver = LQOI.MockLinQuadOptimizer()
        MOIT.lintest(solver, linconfig)
    end

    @testset "Integer Linear tests" begin
        intconfig = MOIT.TestConfig(solve=false)
        solver = LQOI.MockLinQuadOptimizer()
        MOIT.intlineartest(solver, intconfig, ["int2"])
    end

    @testset "ModelLike tests" begin
        intconfig = MOIT.TestConfig()
        solver = LQOI.MockLinQuadOptimizer()
        MOIT.nametest(solver)
        @testset "validtest" begin
            MOIT.validtest(solver)
        end
        @testset "emptytest" begin
            MOIT.emptytest(solver)
        end
        @testset "orderedindicestest" begin
            MOIT.orderedindicestest(solver)
        end
        @testset "canaddconstrainttest" begin
            MOIT.canaddconstrainttest(solver, Float64, Complex{Float64})
        end
        @testset "copytest" begin
            solver2 = LQOI.MockLinQuadOptimizer()
            MOIT.copytest(solver,solver2)
        end
    end
end
