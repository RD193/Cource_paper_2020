import numpy as np
from npplus.fermi import ifd12, fd12
from random import random
import matplotlib.pyplot as plt
from Firs_ascets import *

n = 2e24
NAmND = 1.2e24
p = n + NAmND
x1 = 0.28
x2 = 0.37
Eg1 = 1.424 + 1.247 * x1
Eg2 = 1.424
Eg3 = 1.424 + 1.247 * x2
hi1 = 4.07 - 1.1 * x1
hi2 = 4.07
hi3 = 4.07 - 1.1 * x2
mdn1 = (0.063 + 0.083 * x1) * me
mdn2 = 0.063 * me
mdn3 = (0.063 + 0.083 * x2) * me
mdp1 = (0.51 + 0.25 * x1) * me
mdp2 = 0.51 * me
mdp3 = (0.51 + 0.25 * x2) * me

N = lambda m: 2 * (2 * pi * m * k * T / h ** 2) ** (3 / 2)
Nc = N(mdn1)
Np = N(mdp2)
NP = N(mdp3)
FD1 = NAND1 * np.sqrt(pi) / 2 / Nc
FD2 = NAmND * np.sqrt(pi) / 2 / Np
FD3 = NAND2 * np.sqrt(pi) / 2 / NP


def revFD(v, eta0, eps, bounds=0.1):
    res = fd12(eta0)
    guess = eta0
    while abs(res - v) > eps:
        guess = eta0 + (random() - 0.5) * bounds * 2
        res = fd12(guess)
    return guess


# dFN = 0.7016
# dFp = -0.238
# dFP = -0.376
dFN = ifd12(FD1)
dFp = ifd12(FD2)
dFP = ifd12(FD3)

FN = (-hi1 + dFN * kTeV)
Fp = (-hi2 - Eg2 - dFp * kTeV)
FP = (-hi3 - Eg3 - dFP * kTeV)

F1 = np.ones(11) * FN
F2 = np.ones(11) * Fp
F3 = np.ones(11) * FP

Ec1 = np.ones(11) * (-hi1)
Ec2 = np.ones(11) * (-hi2)
Ec3 = np.ones(11) * (-hi3)

Ev1 = np.ones(11) * (-hi1 - Eg1)
Ev2 = np.ones(11) * (-hi2 - Eg2)
Ev3 = np.ones(11) * (-hi3 - Eg3)

for i in [[xax1, F1, Ec1, Ev1], [xax2, F2, Ec2, Ev2], [xax3, F3, Ec3, Ev3]]:
    plt.plot(i[0], i[1], 'r--')
    plt.plot(i[0], i[2], 'b-')
    plt.plot(i[0], i[3], 'b-')
plt.plot(xax, vaccuum, 'k--')

plt.text(5, 0 - hi1 / 2, 'χ1', color='black', fontsize=8)
plt.text(14, 0 - hi2 / 2, 'χ2', color='black', fontsize=8)
plt.text(23, 0 - hi3 / 2, 'χ3', color='black', fontsize=8)

plt.text(0, 0 - hi1 + 0.2, 'Ec1', color='blue', fontsize=8)
plt.text(11, 0 - hi2 + 0.2, 'Ec2', color='blue', fontsize=8)
plt.text(18, 0 - hi3 + 0.2, 'Ec3', color='blue', fontsize=8)

plt.text(0, FN + 0.1, 'F1', color='red', fontsize=8)
plt.text(11, Fp + 0.1, 'F2', color='red', fontsize=8)
plt.text(22, FP + 0.1, 'F3', color='red', fontsize=8)

plt.text(0, 0 - hi1 - Eg1 + 0.1, 'Ev1', color='blue', fontsize=8)
plt.text(11, 0 - hi2 - Eg2 - 0.2, 'Ev2', color='blue', fontsize=8)
plt.text(21, 0 - hi3 - Eg3 - 0.2, 'Ev3', color='blue', fontsize=8)

plt.text(1.5, -4.6, 'Al(0.28)Ga(0.72)As', color='black', fontsize=8)
plt.text(22.5, -4.6, 'Al(0.37)Ga(0.63)As', color='black', fontsize=8)
plt.text(12.5, -4.8, 'GaAs', color='black', fontsize=8)

plt.text(4.5, -4.9, 'N', color='black', fontsize=9)
plt.text(13.3, -5, 'p', color='black', fontsize=9)
plt.text(24.5, -4.9, 'P', color='black', fontsize=9)

plt.arrow(4.5, 0, 0, -hi1, head_width=0.3, head_length=0.2, linewidth=0.3, color='k', length_includes_head=True)
plt.arrow(4.5, -hi1, 0, hi1, head_width=0.3, head_length=0.2, linewidth=0.3, color='k', length_includes_head=True)

plt.arrow(13.5, 0, 0, -hi2, head_width=0.3, head_length=0.2, linewidth=0.3, color='k', length_includes_head=True)
plt.arrow(13.5, -hi2, 0, hi2, head_width=0.3, head_length=0.2, linewidth=0.3, color='k', length_includes_head=True)

plt.arrow(22.5, 0, 0, -hi3, head_width=0.3, head_length=0.2, linewidth=0.3, color='k', length_includes_head=True)
plt.arrow(22.5, -hi3, 0, hi3, head_width=0.3, head_length=0.2, linewidth=0.3, color='k', length_includes_head=True)

plt.xlabel("x")
plt.ylabel("E, eV")

plt.show()

dEc1 = hi2 - hi1
dEv1 = Eg1 - Eg2 - dEc1

eps1 = (13.1 - 3.0 * x1) * eps0
eps2 = 13.1 * eps0
K1 = 1 + eps2 * NAmND / (eps1 * NAmND)

VD1 = FN - Fp
VDp1 = VD1 / K1
VDN = VD1 - VDp1

xN = np.sqrt(2 * eps1 * VDN / (e * NAND1))
xp1 = np.sqrt(2 * eps2 * VDp1 / (e * NAmND))

xaxN = np.linspace(-xN, 0, 20)
xaxp1 = np.linspace(0, xp1, 20)
EvNp = np.append(VDN - e * NAND1 / (2 * eps1) * (xN ** 2 - (xN + xaxN) ** 2),
                 VD1 + dEv1 - (e * NAmND / (2 * eps2)) * (xp1 - xaxp1) ** 2)
EcNp = EvNp + np.append(Eg1 * np.ones(20), Eg2 * np.ones(20))

dEc2 = hi2 - hi3
dEv2 = Eg3 - Eg2 - dEc2

eps3 = (13.1 - 3.0 * x2) * eps0
K2 = 1 + eps3 * NAmND / (eps2 * NAND2)

VD2 = Fp - FP
VDp2 = VD2 / K2
VDP = VD2 - VDp2

xp2 = np.sqrt(2 * eps2 * VDp2 / (e * NAmND))
xP = np.sqrt(2 * eps3 * VDP / (e * NAND2))

xaxp2 = np.linspace(-xp2, 0, 20)
xaxP = np.linspace(0, xP, 20)
EvpP = np.append(VDp2 - e * NAmND / (2 * eps2) * (xp2 ** 2 - (xp2 + xaxp2) ** 2),
                 VDp2 - dEv2 + (e * NAND2 / (2 * eps3)) * (xP ** 2 - (xP - xaxP) ** 2)) + VD1 + dEv1
EcpP = EvpP + np.append(Eg2 * np.ones(20), Eg3 * np.ones(20))

plt.figure(2)
plt.axvline(0, 0, Eg3 + VDp2 - dEv2 + (e * NAND2 / (2 * eps3)) * (xP ** 2) + VD1 + dEv1, color='k', linestyle='--',
            linewidth=0.5)
plt.axvline(d, 0, Eg3 + VDp2 - dEv2 + (e * NAND2 / (2 * eps3)) * (xP ** 2) + VD1 + dEv1, color='k', linestyle='--',
            linewidth=0.5)

plt.plot(np.append(xaxN, xaxp1), EvNp, 'b-')
plt.plot(np.append(xaxN, xaxp1), EcNp, 'b-')
plt.plot(np.append(xaxp2 + d, xaxP + d), EvpP, 'b-')
plt.plot(np.append(xaxp2 + d, xaxP + d), EcpP, 'b-')

plt.plot([-xN - 1e-7, xP + d + 1e-7], np.ones(2) * (Eg1 + dFN * kTeV), 'g--')
plt.text(-xN - 1e-7, Eg1 - dFN * kTeV - 0.25, 'F', color='g', fontsize=10)

plt.plot([xp1, d - xp2], (VD1 + dEv1) * np.ones(2), 'b-')
plt.plot([xp1, d - xp2], (VD1 + dEv1 + Eg2) * np.ones(2), 'b-')
plt.plot([-2e-7, -xN], np.zeros(2), 'b-')
plt.text(-1.8e-7, 0.1, 'Ev', color='b', fontsize=10)
plt.plot([-2e-7, -xN], np.zeros(2) + Eg1, 'b-')
plt.text(-1.8e-7, Eg1 + 0.1, 'Ec', color='b', fontsize=10)
plt.plot([d + xP, d + xP + 1e-7], np.zeros(2) + VDp2 - dEv2 + (e * NAND2 / (2 * eps3)) * (xP ** 2) + VD1 + dEv1, 'b-')
plt.plot([d + xP, d + xP + 1e-7], np.zeros(2) + Eg3 + VDp2 - dEv2 + (e * NAND2 / (2 * eps3)) * (xP ** 2) + VD1 + dEv1,
         'b-')

plt.text(-1e-7, 3.7, 'Al0.28Ga072As', color='k', fontsize=10)
plt.text(5e-8, 3.7, 'GaAs', color='k', fontsize=10)
plt.text(d + 0.11e-7, 3.7, 'Al0.37Ga063As', color='k', fontsize=10)
plt.xlabel('x')
plt.ylabel('E, eV')
plt.show()

NN2 = N(mdn3)
NP2 = N(mdp1)
dFN1 = ifd12(NAND1 * np.sqrt(pi) / 2 / Nc) * kTeV
dFc_ = ifd12(n * np.sqrt(pi) / 2 / NN2) * kTeV
dFP1 = ifd12(NAND2 * np.sqrt(pi) / 2 / NP) * kTeV
dFv_ = ifd12((NAmND + n) * np.sqrt(pi) / 2 / NP2) * kTeV

Fc_ = -hi2 + dFc_
Fv_ = -hi2 - Eg2 - dFv_

volt = -Fv_ + Fc_

VDp_VD1 = (VD1 - volt) / K1
VDN_VD1 = (VD1 - volt) * (1 - 1 / K1)

xN = np.sqrt((2 * eps1 * VDN_VD1) / (e * NAND1))
xp1 = np.sqrt((2 * eps2 * VDp_VD1) / (e * NAmND))
xaxN = np.linspace(-xN, 0, 20)
xaxp1 = np.linspace(0, xp1, 20)

EvNp = np.append(VDN_VD1 - e * NAND1 / (2 * eps1) * (xN ** 2 - (xN + xaxN) ** 2),
                 VD1 - volt + dEv1 - (e * NAmND / (2 * eps2)) * (xp1 - xaxp1) ** 2)
EcNp = EvNp + np.append(Eg1 * np.ones(20), Eg2 * np.ones(20))


xp2 = np.sqrt(2 * eps2 * VDp2 / (e * NAmND))
xP = np.sqrt(2 * eps3 * VDP / (e * NAND2))
xaxp2 = np.linspace(-xp2, 0, 20)
xaxP = np.linspace(0, xP, 20)

EvpP = np.append(VD2/K2- volt - e * NAmND / (2 * eps2) * (xp2 ** 2 - (xp2 + xaxp2) ** 2),
                 VD2/K2 -dEv2 - volt + (e * NAND2 / (2 * eps3)) * (xP ** 2 - (xP - xaxP) ** 2)) + VD1 + dEv1
EcpP = EvpP + np.append(Eg2 * np.ones(20), Eg3 * np.ones(20))

plt.figure(3)
plt.axvline(0, 0, Eg3 + VDp2 - dEv2 + (e * NAND2 / (2 * eps3)) * (xP ** 2) + VD1 + dEv1, color='k', linestyle='--',
            linewidth=0.5)
plt.axvline(d, 0, Eg3 + VDp2 - dEv2 + (e * NAND2 / (2 * eps3)) * (xP ** 2) + VD1 + dEv1, color='k', linestyle='--',
            linewidth=0.5)

plt.plot(np.append(xaxN, xaxp1), EvNp, 'b-')
plt.plot(np.append(xaxN, xaxp1), EcNp, 'b-')
plt.plot(np.append(xaxp2 + d, xaxP + d), EvpP, 'b-')
plt.plot(np.append(xaxp2 + d, xaxP + d), EcpP, 'b-')
plt.plot([-1.2e-7, d], np.ones(2) * (Eg1 + dFN*kTeV), 'g--')
plt.plot([0, d + xP + 1e-7], np.ones(2) * (Eg1 + dFN*kTeV - volt), 'g--')

plt.text(-xN - 2e-8, Eg1 + dFN * kTeV - 0.15, 'Fc*', color='g', fontsize=10)
plt.text(d + xP + 2e-8, Eg1 + dFN*kTeV - volt + 0.1, 'Fv*', color='g', fontsize=10)

plt.plot([xp1, d - xp2], (VD1 + dEv1) * np.ones(2) - volt, 'b-')
plt.plot([xp1, d - xp2], (VD1 + dEv1 + Eg2) * np.ones(2) - volt, 'b-')
plt.plot([-1.2e-7, -xN], np.zeros(2), 'b-')
plt.text(-1e-7, 0.1, 'Ev', color='b', fontsize=10)
plt.plot([-1.2e-7, -xN], np.zeros(2) + Eg1, 'b-')
plt.text(-1e-7, Eg1 + 0.1, 'Ec', color='b', fontsize=10)
plt.plot([d + xP, d + xP + 1e-7], np.zeros(2) + VDp2 - dEv2 + (e * NAND2 / (2 * eps3)) * (xP ** 2) + VD1 + dEv1 - volt, 'b-')
plt.plot([d + xP, d + xP + 1e-7], np.zeros(2) + Eg3 + VDp2 - dEv2 + (e * NAND2 / (2 * eps3)) * (xP ** 2) + VD1 + dEv1 - volt,
         'b-')

plt.text(-1e-7, 2.4, 'Al0.28Ga072As', color='k', fontsize=10)
plt.text(3e-8, 2.4, 'GaAs', color='k', fontsize=10)
plt.text(d + 0.11e-7, 2.4, 'Al0.37Ga063As', color='k', fontsize=10)