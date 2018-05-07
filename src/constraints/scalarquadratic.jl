#=
        ScalarQuadraticFunction -in- LessThan
        ScalarQuadraticFunction -in- GreaterThan
        ScalarQuadraticFunction -in- EqualTo

    TODO:
        - constraint sets
        - constraint functions
        - deleting constraints
=#
constrdict(m::LinQuadOptimizer, ::QCI{LE})  = cmap(m).q_less_than
constrdict(m::LinQuadOptimizer, ::QCI{GE})  = cmap(m).q_greater_than
constrdict(m::LinQuadOptimizer, ::QCI{EQ})  = cmap(m).q_equal_to

function MOI.addconstraint!(m::LinQuadOptimizer, func::Quad, set::S) where S <: Union{LE, GE, EQ}
    addquadraticconstraint!(m, func, set)
    m.last_constraint_reference += 1
    ref = QCI{S}(m.last_constraint_reference)
    dict = constrdict(m, ref)
    push!(m.qconstraint_primal_solution, NaN)
    push!(m.qconstraint_dual_solution, NaN)
    # dict[ref] = get_number_quadratic_constraints(m)
    dict[ref] = length(m.qconstraint_primal_solution)
    return ref
end

function addquadraticconstraint!(m::LinQuadOptimizer, func::Quad, set::S) where S<: Union{LE, GE, EQ}
    addquadraticconstraint!(m, func, lqs_char(m,set), _getrhs(set))
end

function addquadraticconstraint!(m::LinQuadOptimizer, f::Quad, sense::Cchar, rhs::Float64)
    if abs(f.constant) > 0
        warn("Constant in quadratic function. Moving into set")
    end
    ri, ci, vi = reduceduplicates(
        getcol.(m, f.quadratic_rowvariables),
        getcol.(m, f.quadratic_colvariables),
        f.quadratic_coefficients
    )
    add_quadratic_constraint!(m,
        getcol.(m, f.affine_variables),
        f.affine_coefficients,
        rhs - f.constant,
        sense,
        ri, ci, vi
    )
end

function reduceduplicates(rowi::Vector{T}, coli::Vector{T}, vals::Vector{S}) where T where S
    @assert length(rowi) == length(coli) == length(vals)
    d = Dict{Tuple{T, T},S}()
    for (r,c,v) in zip(rowi, coli, vals)
        if haskey(d, (r,c))
            d[(r,c)] += v
        else
            d[(r,c)] = v
        end
    end
    ri = Vector{T}(length(d))
    ci = Vector{T}(length(d))
    vi = Vector{S}(length(d))
    for (i, (key, val)) in enumerate(d)
        ri[i] = key[1]
        ci[i] = key[2]
        vi[i] = val
    end
    ri, ci, vi
end
