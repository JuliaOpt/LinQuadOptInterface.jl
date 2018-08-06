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
