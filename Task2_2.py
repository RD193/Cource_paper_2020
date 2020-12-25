import math
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from sympy.abc import k

da = 6e-9
q0 = 1.6e-19
me = 9.11e-31
mn = 0.067*me
mp = 0.48*me
ev = 6.242e18
V0 = 55e-3
NAmND = 1.2e24
n = 2e24
h =1.0510e-34
kT = 4.11e-21

x = 0.158
Eg1 = 1.424 + 1.247*x
Eg2 = 1.424


n2d = n*da

Efn = kT*np.log(np.exp(n2d*np.pi*h**2/(mn*kT)) - 1)*ev
print(Efn)
p2d = (NAmND + n)*da

Efp = kT*np.log(np.exp(p2d*np.pi*h**2/(mp*kT)) - 1)*ev
print(Efp)

hi1 = 4.07 - 1.1*x
hi2 = 4.07

ye = 25.817057012959342*pow(10,-3)
yh = 14.898742709110866*pow(10,-3)
mainx = [-9, -3, -3, 3, 3, 9]
mainyc = [-hi1, -hi1, -hi1-V0, -hi1-V0, -hi1, -hi1]
mainyv = [-hi1-Eg1, -hi1-Eg1, -hi1-V0-Eg2, -hi1-V0-Eg2, -hi1-Eg1, -hi1-Eg1]
plt.plot(mainx, mainyc,'k')
plt.plot(mainx, mainyv,'k')

plt.plot(mainx[2:4], [-hi1-V0 +Efn, -hi1-V0 +Efn],'r--', linewidth =0.7)
plt.plot(mainx[2:4], [-hi1-V0-Eg2 + Efp, -hi1-V0-Eg2+ Efp ], 'r--', linewidth =0.7)

plt.plot(mainx[2:4], [-hi1-V0 +ye, -hi1-V0 +ye], 'b--', linewidth =0.7)
plt.plot(mainx[2:4], [-hi1-V0-Eg2 -yh, -hi1-V0-Eg2-yh ], 'b--', linewidth =0.7)




plt.text(-2,-hi1-V0 +ye-0.1, 'En', color='blue', fontsize=8)
plt.text(-2,-hi1-V0-Eg2 -yh+0.05, 'Eh', color='blue', fontsize=8)
plt.text(2.05,-hi1-V0-Eg2 + Efp+ 0.05, 'Fv*', color='red', fontsize=8)
plt.text(2.05,-hi1-V0 +Efn-0.09, 'Fc*', color='red', fontsize=8)

plt.xlabel("x, nm")
plt.show()


















