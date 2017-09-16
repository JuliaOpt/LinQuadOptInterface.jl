# Main

# LinQuadSolver # Abstract type
function lqs_setparam!(env, name, val) end
function lqs_setlogfile!(env, path) end
function lqs_getprobtype(m) end

# Constraints

function lqs_chgbds!(m, colvec, valvec, sensevec) end
function lqs_getlb(m, col) end
function lqs_getub(m, col) end
function lqs_getnumrows(m) end
function lqs_addrows!(m, rowvec, colvec, coefvec, sensevec, rhsvec) end
function lqs_getrhs(m, row) end
function lqs_getrows(m, rowvec) end
function lqs_getcoef(m, row, col) end #??
function lqs_chgcoef!(m, row, col, coef) end
function lqs_delrows!(m, row, row2) end
function lqs_chgctype!(m, colvec, typevec) end
function lqs_chgsense!(m, rowvec, sensevec) end

function lqs_vartype_map(m) end
function lqs_make_problem_type_integer(m) end
function lqs_make_problem_type_continuous(m) end

function lqs_addsos!(m, colvec, valvec, typ) end
function lqs_delsos!(m, idx, idx2) end
function lqs_sertype_map(m) end
function lqs_getsos(m, idx) end

function lqs_getnumqconstrs(m) end
function lqs_addqconstr!(m, cols,coefs,rhs,sense, I,J,V) end

function lqs_chgrngval!(m, rows, vals) end# later
function lqs_ctrtype_map(m) end

#Objective

function lqs_copyquad!(m, intvec,intvec2, floatvec) end#?
function lqs_chgobj!(m, colvec,coefvec) end
function lqs_chgobjsen!(m, symbol) end
function lqs_getobj(m) end
function lqs_getobjsen(m) end

#Solve

function lqs_mipopt!(m) end
function lqs_qpopt!(m) end
function lqs_lpopt!(m) end
function lqs_getstat(m) end
function lqs_solninfo(m) end # complex
function lqs_getx!(m, place) end
function lqs_getax!(m, place) end
function lqs_getdj!(m, place) end
function lqs_getpi!(m, place) end

function lqs_getobjval(m) end
function lqs_getbestobjval(m) end
function lqs_getmiprelgap(m) end
function lqs_getitcnt(m) end
function lqs_getbaritcnt(m) end
function lqs_getnodecnt(m) end

function lqs_dualfarkas!(m, place) end
function lqs_getray!(m, place) end

function lqs_terminationstatus end
function lqs_primalstatus end
function lqs_dualstatus end

# Variables

function lqs_getnumcols(m) end
function lqs_newcols!(m, int) end
function lqs_delcols!(m, col, col2) end
function lqs_addmipstarts!(m, colvec, valvec) end