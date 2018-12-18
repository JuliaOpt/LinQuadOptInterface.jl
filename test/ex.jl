using Xpress, Base.Test, MathOptInterface, MathOptInterface.Test, MathOptInterfaceXpress
using CPLEX, MathOptInterfaceCPLEX
using Gurobi, MathOptInterfaceGurobi
using GLPK, MathOptInterfaceGLPK


const MOI = MathOptInterface

const MOIU = MathOptInterface.Utilities

solvers = [XpressOptimizer(),
          CPLEXOptimizer(),
          GLPKOptimizerLP(),
        #   GLPKOptimizerMIP(),
          GurobiOptimizer()]


for solver in solvers
    m = Model()
    x = @variable(m)
    y = @variable(m)
    @objective(m, Min, -x)

    c = @constraint(m, x + y <= 1)
    c1 = @constraint(m, y >= 0)
    c2 = @constraint(m, x >= 0)

    MOIU.resetoptimizer!(m, solver)
    MOIU.attachoptimizer!(m)

    JuMP.optimize(m)

    @test JuMP.hasresultvalues(m)

    @test JuMP.terminationstatus(m) == MOI.Success
    @test JuMP.primalstatus(m) == MOI.FEASIBLE_POINT

    JuMP.resultvalue(x)
    JuMP.resultvalue(y)
    JuMP.resultvalue(x + y)
    JuMP.objectivevalue(m)

    JuMP.dualstatus(m)
    JuMP.resultdual(c)
end