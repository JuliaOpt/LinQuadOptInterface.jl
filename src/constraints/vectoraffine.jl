#=
    Vector valued constraints
=#
constrdict(m::LinQuadOptimizer, ::VLCI{MOI.Nonnegatives}) = cmap(m).nonnegatives
constrdict(m::LinQuadOptimizer, ::VLCI{MOI.Nonpositives}) = cmap(m).nonpositives
constrdict(m::LinQuadOptimizer, ::VLCI{MOI.Zeros})        = cmap(m).zeros

function MOI.addconstraint!(m::LinQuadOptimizer, func::VecLin, set::S) where S <: VecLinSets
    @assert MOI.dimension(set) == length(func.constant)

    nrows = get_number_linear_constraints(m)
    addlinearconstraint!(m, func, lqs_char(m,set))
    nrows2 = get_number_linear_constraints(m)

    m.last_constraint_reference += 1
    ref = VLCI{S}(m.last_constraint_reference)

    dict = constrdict(m, ref)
    dict[ref] = collect(nrows+1:nrows2)
    for i in 1:MOI.dimension(set)
        push!(m.constraint_primal_solution, NaN)
        push!(m.constraint_dual_solution, NaN)
        push!(m.constraint_constant, func.constant[i])
    end
    ref
end

function addlinearconstraint!(m::LinQuadOptimizer, func::VecLin, sense::Cchar)
    @assert length(func.outputindex) == length(func.variables) == length(func.coefficients)
    # get list of unique rows
    rows = unique(func.outputindex)
    @assert length(rows) == length(func.constant)
    # sort into row order
    pidx = sortperm(func.outputindex)
    cols = getcol.(m, func.variables)[pidx]
    vals = func.coefficients[pidx]
    # loop through to gte starting position of each row
    rowbegins = Vector{Int}(length(rows))
    rowbegins[1] = 1
    cnt = 1
    for i in 2:length(pidx)
        if func.outputindex[pidx[i]] != func.outputindex[pidx[i-1]]
            cnt += 1
            rowbegins[cnt] = i
        end
    end
    add_linear_constraints!(m, rowbegins, cols, vals, fill(sense, length(rows)), -func.constant)
end

MOI.canmodifyconstraint(m::LinQuadOptimizer, ::VLCI{<: VecLinSets}, ::Type{MOI.VectorConstantChange{Float64}}) = true
function MOI.modifyconstraint!(m::LinQuadOptimizer, ref::VLCI{<: VecLinSets}, chg::MOI.VectorConstantChange{Float64})
    @assert length(chg.new_constant) == length(m[ref])
    for (r, v) in zip(m[ref], chg.new_constant)
        change_coefficient!(m, r, 0, -v)
        m.constraint_constant[r] = v
    end
end

#=
    Constraint set of Linear function
=#

MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintSet, ::Type{VLCI{S}}) where S <: VecLinSets = true
function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintSet, c::VLCI{S}) where S <: VecLinSets
    constraint_indices = m[c]
    S(length(constraint_indices))
end

#=
    Constraint function of Linear function
=#

MOI.canget(m::LinQuadOptimizer, ::MOI.ConstraintFunction, ::Type{VLCI{S}}) where S <: VecLinSets = true
function MOI.get(m::LinQuadOptimizer, ::MOI.ConstraintFunction, c::VLCI{<: VecLinSets})
    ctrs = m[c]
    n = length(ctrs)
    out = MOI.VectorAffineFunction(Int[],VarInd[],Float64[],Float64[])
    for i in 1:n
        rhs = get_rhs(m, ctrs[i])
        push!(out.constant, -rhs)

        # TODO more efficiently
        colidx, coefs = get_linear_constraint(m, ctrs[i])
        append!(out.variables, m.variable_references[colidx+1])
        append!(out.coefficients, coefs)
        append!(out.outputindex, i*ones(Int,length(coefs)))
    end
    return out
end
