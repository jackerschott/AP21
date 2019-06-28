import datproc.print as dpr

import stokes as st
import hagen_poiseuille as hp

print(dpr.val(st.eta, st.d_eta, name='η_st', unit='Pa s'))
print(dpr.val(hp.eta, hp.d_eta, name='η_hp', unit='Pa s'))
print(dpr.dev(st.eta, st.d_eta, hp.eta, hp.d_eta, name='η_st, η_hp'))

Re_crit = 3
d_Re_crit = 1
Re_crit_theo = 1
print(dpr.val(Re_crit, d_Re_crit, name='Re_crit'))
print(dpr.val(Re_crit_theo, name='Re_crit_theo'))
print(dpr.dev(Re_crit, d_Re_crit, Re_crit_theo, name='Re_crit, Re_crit_theo'))
