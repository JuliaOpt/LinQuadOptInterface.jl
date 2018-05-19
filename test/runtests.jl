#push!(Base.LOAD_PATH,joinpath(dirname(@__FILE__),"..",".."))

using Base.Test, MathOptInterface
using LinQuadOptInterface

const MOIT = MathOptInterface.Test
const LQOI = LinQuadOptInterface


@testset "Linear tests" begin
    linconfig = MOIT.TestConfig(solve = false)
    solver = LQOI.MockLinQuadOptimizer()
    MOIT.contlineartest(solver , linconfig, ["linear10"])
end