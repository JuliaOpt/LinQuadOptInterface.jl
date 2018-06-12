#push!(Base.LOAD_PATH,joinpath(dirname(@__FILE__),"..",".."))

using Base.Test, MathOptInterface
using LinQuadOptInterface

const MOIT = MathOptInterface.Test
const LQOI = LinQuadOptInterface


@testset "LinQuadOptInterface" begin
    @testset "Linear tests" begin
        linconfig = MOIT.TestConfig(solve = false)
        solver = LQOI.MockLinQuadOptimizer()
        MOIT.contlineartest(solver , linconfig, ["linear10"])
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
        intconfig = MOIT.TestConfig(solve=false)
        solver = LQOI.MockLinQuadOptimizer()
        MOIT.validtest(solver)
        MOIT.emptytest(solver)
        solver2 = LQOI.MockLinQuadOptimizer()
        MOIT.copytest(solver,solver2)
    end
end
;