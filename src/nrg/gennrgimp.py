#!/usr/bin/env python

# This script generates nrgimp.f90.

alg = ["standard", "celldec"]   # algos to calculate energies. nrg_algo.f90 must exist.
box = ["generic", "cubic", "cubicXY", "hexagonal"] # box types. function 'dist2_mi_box' must exist.

deftmpl = """
#define DIST2_MI dist2_mi_{1}
#define NRG_1P  nrg_1p_{0}_{1}
#define NRG_BOX  nrg_box_{0}_{1}
#include "nrg_{0}.f90"
#undef DIST2_MI
#undef NRG_1P
#undef NRG_BOX
"""

defs = ""
privs_nrg = ""
privs_box = "private :: pbc,"
nrg_sel = open("nrginterf.f90").read()

cases_1p = ""
cases_box = ""
cases_tmpl = "    case({2})\n      nrg_{3} = nrg_{3}_{0}_{1}({4})\n"
magn = 100 # or len(alg)

for j in range(0,len(box)):
    privs_box += "dist2_mi_{0}".format(box[j])
    if j < len(box)-1: privs_box += ", "

for i in range(0, len(alg)):
    a = alg[i]
    for j in range(0,len(box)):
        b = box[j]
        defs += deftmpl.format(a,b)
        privs_nrg += "private :: nrg_1p_{0}_{1}, nrg_box_{0}_{1}\n".format(a,b)
        ic = magn * i + j
        cases_1p += cases_tmpl.format(a, b, ic, "1p", "st,b,par")
        cases_box += cases_tmpl.format(a, b, ic, "box", "st,b,press")

with open("nrgprivs.f90", "w") as f:
    f.write(privs_box)
    f.write("\n")
    f.write(privs_nrg)
    f.write("\n")
with open("nrgimp.f90", "w") as f:
    f.write("include \"dist2_mi.f90\"\n\n")
    f.write(nrg_sel.format(cases_1p, cases_box))
    f.write(defs)


