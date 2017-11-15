# Main

# LinQuadSolver # Abstract type
function lqs_setparam!(env, name, val) end
function lqs_setlogfile!(env, path) end
function lqs_getprobtype(m::LinQuadSolverInstance) end

lqs_supported_constraints(s) = []
lqs_supported_objectives(s) = []

# Constraints

function lqs_chgbds!(m::LinQuadSolverInstance, colvec, valvec, sensevec) end
function lqs_getlb(m::LinQuadSolverInstance, col) end
function lqs_getub(m::LinQuadSolverInstance, col) end
function lqs_getnumrows(m::LinQuadSolverInstance) end
function lqs_addrows!(m::LinQuadSolverInstance, rowvec, colvec, coefvec, sensevec, rhsvec) end
function lqs_getrhs(m::LinQuadSolverInstance, row) end
function lqs_getrows(m::LinQuadSolverInstance, rowvec) end
function lqs_getcoef(m::LinQuadSolverInstance, row, col) end #??
function lqs_chgcoef!(m::LinQuadSolverInstance, row, col, coef) end
function lqs_delrows!(m::LinQuadSolverInstance, row, row2) end
function lqs_chgctype!(m::LinQuadSolverInstance, colvec, typevec) end
function lqs_chgsense!(m::LinQuadSolverInstance, rowvec, sensevec) end

function lqs_vartype_map(m::LinQuadSolverInstance) end
function lqs_make_problem_type_integer(m::LinQuadSolverInstance) end
function lqs_make_problem_type_continuous(m::LinQuadSolverInstance) end

function lqs_addsos!(m::LinQuadSolverInstance, colvec, valvec, typ) end
function lqs_delsos!(m::LinQuadSolverInstance, idx, idx2) end
function lqs_sertype_map(m::LinQuadSolverInstance) end
function lqs_getsos(m::LinQuadSolverInstance, idx) end

function lqs_getnumqconstrs(m::LinQuadSolverInstance) end
function lqs_addqconstr!(m::LinQuadSolverInstance, cols,coefs,rhs,sense, I,J,V) end

function lqs_chgrngval!(m::LinQuadSolverInstance, rows, vals) end# later
function lqs_ctrtype_map(m::LinQuadSolverInstance) end

#Objective

function lqs_copyquad!(m::LinQuadSolverInstance, intvec,intvec2, floatvec) end#?
function lqs_chgobj!(m::LinQuadSolverInstance, colvec,coefvec) end
function lqs_chgobjsen!(m::LinQuadSolverInstance, symbol) end
function lqs_getobj(m::LinQuadSolverInstance) end
function lqs_getobjsen(m::LinQuadSolverInstance) end

#Solve

function lqs_mipopt!(m::LinQuadSolverInstance) end
function lqs_qpopt!(m::LinQuadSolverInstance) end
function lqs_lpopt!(m::LinQuadSolverInstance) end
function lqs_getstat(m::LinQuadSolverInstance) end
function lqs_solninfo(m::LinQuadSolverInstance) end # complex
function lqs_getx!(m::LinQuadSolverInstance, place) end
function lqs_getax!(m::LinQuadSolverInstance, place) end
function lqs_getdj!(m::LinQuadSolverInstance, place) end
function lqs_getpi!(m::LinQuadSolverInstance, place) end

function lqs_getobjval(m::LinQuadSolverInstance) end
function lqs_getbestobjval(m::LinQuadSolverInstance) end
function lqs_getmiprelgap(m::LinQuadSolverInstance) end
function lqs_getitcnt(m::LinQuadSolverInstance) end
function lqs_getbaritcnt(m::LinQuadSolverInstance) end
function lqs_getnodecnt(m::LinQuadSolverInstance) end

function lqs_dualfarkas!(m::LinQuadSolverInstance, place) end
function lqs_getray!(m::LinQuadSolverInstance, place) end

function lqs_terminationstatus end
function lqs_primalstatus end
function lqs_dualstatus end

# Variables

function lqs_getnumcols(m::LinQuadSolverInstance) end
function lqs_newcols!(m::LinQuadSolverInstance, int) end
function lqs_delcols!(m::LinQuadSolverInstance, col, col2) end
function lqs_addmipstarts!(m::LinQuadSolverInstance, colvec, valvec) end