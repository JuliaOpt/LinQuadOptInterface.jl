Main

LinQuadSolver # Abstract type
lqs_setparam!(env, name, val)
lqs_setlogfile!(env, path)
lqs_getprobtype(m)

Constraints

lqs_chgbds!(m, colvec, valvec, sensevec)
lqs_getlb(m, col)
lqs_getub(m, col)
lqs_getnumrows(m)
lqs_addrows!(m, rowvec, colvec, coefvec, sensevec, rhsvec)
lqs_getrhs(m, row)
colvec, coef = lqs_getrows(m, rowvec)
lqs_getcoef(m, row, col) #??
lqs_chgcoef!(m, row, col, coef)
lqs_delrows!(m, row, row)
lqs_chgctype!(m, colvec, typevec)
lqs_chgsense!(m, rowvec, sensevec)

# lqs_addsos(m, colvec, valvec, typ)
# lqs_delsos(m, idx, idx)
# lqs_getsos(m, idx)

# lqs_getnumqconstrs(m)
# lqs_addqconstr(m, cols,coefs,rhs,sense, I,J,V)

# lqs_chgrngval # later

Objective

lqs_copyquad(m, intvec,intvec, floatvec) #?
lqs_chgobj!(m, colvec,coefvec)
lqs_chgobjsen(!m, symbol)
lqs_getobj(m)

Solve

lqs_mipopt!(m)
lqs_qpopt!(m)
lqs_lpopt!(m)
lqs_getstat(m)
lqs_solninfo(m) # complex
lqs_getx!(m, place)
lqs_getax!(m, place)
lqs_getdj!(m, place)
lqs_getpi!(m, place)

lqs_getobjval(m)
lqs_getbestobjval(m)
lqs_getmiprelgap(m)
lqs_getitcnt(m)
lqs_getbaritcnt(m)
lqs_getnodecnt(m)

lqs_termination_status_map(m) # = TERMINATION_STATUS_MAP
lqs_sol_basic(m) #
lqs_sol_nonbasic(m)
lqs_sol_primal(m)
lqs_sol_none(m)

# lqs_dualopt(m)
# lqs_dualfarkas(m, place)
# lqs_getray(m, place)

Variables

lqs_getnumcols(m)
lqs_newcols!(m, int)
lqs_delcols!(m, col, col)
# lqs_addmipstarts(m, colvec, valvec)