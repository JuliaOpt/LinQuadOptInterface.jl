# Main

"""
    LinQuadModel(M,env)

Initializes a model given a model type M and a env, that might be a nothing for some solvers.
"""
function LinQuadModel end

# LinQuadSolver # Abstract type
function lqs_setparam!(env, name, val) end
function lqs_setlogfile!(env, path) end
function lqs_getprobtype(m::LinQuadOptimizer) end

lqs_supported_constraints(s) = []
lqs_supported_objectives(s) = []

# Constraints

function lqs_chgbds!(m::LinQuadOptimizer, colvec, valvec, sensevec) end
function lqs_getlb(m::LinQuadOptimizer, col) end
function lqs_getub(m::LinQuadOptimizer, col) end
function lqs_getnumrows(m::LinQuadOptimizer) end
function lqs_addrows!(m::LinQuadOptimizer, rowvec, colvec, coefvec, sensevec, rhsvec) end
function lqs_getrhs(m::LinQuadOptimizer, row) end
function lqs_getrows(m::LinQuadOptimizer, rowvec) end
function lqs_getcoef(m::LinQuadOptimizer, row, col) end #??
function lqs_chgcoef!(m::LinQuadOptimizer, row, col, coef) end
function lqs_delrows!(m::LinQuadOptimizer, row, row2) end
function lqs_chgctype!(m::LinQuadOptimizer, colvec, typevec) end
function lqs_chgsense!(m::LinQuadOptimizer, rowvec, sensevec) end

function lqs_vartype_map(m::LinQuadOptimizer) end
function lqs_make_problem_type_integer(m::LinQuadOptimizer) end
function lqs_make_problem_type_continuous(m::LinQuadOptimizer) end

function lqs_addsos!(m::LinQuadOptimizer, colvec, valvec, typ) end
function lqs_delsos!(m::LinQuadOptimizer, idx, idx2) end
function lqs_sertype_map(m::LinQuadOptimizer) end
function lqs_getsos(m::LinQuadOptimizer, idx) end

function lqs_getnumqconstrs(m::LinQuadOptimizer) end
function lqs_addqconstr!(m::LinQuadOptimizer, cols,coefs,rhs,sense, I,J,V) end

function lqs_chgrngval!(m::LinQuadOptimizer, rows, vals) end# later
function lqs_ctrtype_map(m::LinQuadOptimizer) end

#Objective

function lqs_copyquad!(m::LinQuadOptimizer, intvec,intvec2, floatvec) end#?
function lqs_chgobj!(m::LinQuadOptimizer, colvec,coefvec) end
function lqs_chgobjsen!(m::LinQuadOptimizer, symbol) end
function lqs_getobj(m::LinQuadOptimizer) end
function lqs_getobjsen(m::LinQuadOptimizer) end

#Solve

function lqs_mipopt!(m::LinQuadOptimizer) end
function lqs_qpopt!(m::LinQuadOptimizer) end
function lqs_lpopt!(m::LinQuadOptimizer) end
function lqs_getstat(m::LinQuadOptimizer) end
function lqs_solninfo(m::LinQuadOptimizer) end # complex
function lqs_getx!(m::LinQuadOptimizer, place) end
function lqs_getax!(m::LinQuadOptimizer, place) end
function lqs_getqcax!(m::LinQuadOptimizer, place) end
function lqs_getdj!(m::LinQuadOptimizer, place) end
function lqs_getpi!(m::LinQuadOptimizer, place) end
function lqs_getqcpi!(m::LinQuadOptimizer, place) end

function lqs_getobjval(m::LinQuadOptimizer) end
function lqs_getbestobjval(m::LinQuadOptimizer) end
function lqs_getmiprelgap(m::LinQuadOptimizer) end
function lqs_getitcnt(m::LinQuadOptimizer) end
function lqs_getbaritcnt(m::LinQuadOptimizer) end
function lqs_getnodecnt(m::LinQuadOptimizer) end

function lqs_dualfarkas!(m::LinQuadOptimizer, place) end
function lqs_getray!(m::LinQuadOptimizer, place) end

function lqs_terminationstatus end
function lqs_primalstatus end
function lqs_dualstatus end

# Variables

function lqs_getnumcols(m::LinQuadOptimizer) end
function lqs_newcols!(m::LinQuadOptimizer, int) end
function lqs_delcols!(m::LinQuadOptimizer, col, col2) end
function lqs_addmipstarts!(m::LinQuadOptimizer, colvec, valvec) end
