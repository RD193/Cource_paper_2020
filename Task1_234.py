from numpy import sqrt, log, pi, arcsin, linspace
import matplotlib.pyplot as plt

h = 4.135e-15
hSI = 6.626e-34
c = 3e8
e = 1.6e-19
k = 1.38e-23
T = 298
lam = 884e-9
aomg = 22930
bomg = 1e3
B = 7.1e-16
ksi = 1
L = 1.105e-3
W = 12e-6
n0 = 5e23
n = 2e24
da = 1.1e-6
E = h*c /lam
E0 = 1.435
E1 = 2.9
E2 = 5
A = 0.584
G1 = 20.0432
G2 = 151.197

A0 = aomg/(n - n0)

m = sqrt(1 + A/pi * log((E1 **2 - E**2)/(E0**2 - E**2)) + G1/(E1**2 - E**2) + G2/(E2**2 - E**2))
d0 = lam / 2/m

R = ((m -1) / (m + 1))**2
a0 = A0*n0
Lmin = log(1/R)/(aomg - bomg)


J = e * B *da / (A0 * ksi)**2 *(a0*ksi + bomg + 1/L*log(1/R))**2
I = J*W*L

q = 2*m*L/lam
dlamq = lam**2 / (2 * m * L)
dnuq = c / (2 * m * L)
Q = 4*L*m*pi/(lam*(1 - R))
FWHM = lam / Q

dnuspont = 3*k*T/(2*hSI)
dnulas = ((dnuspont**2)*dnuq)**(1/3)

N = dnulas/dnuq

Thetall = lam / da
Theta_l_ = 2 * arcsin(lam / W)

eta_opt = log(R) / (log(R) - L*bomg)

xax = linspace(5e7,5e9, 100)
eta = eta_opt*( 1 - J/xax)
plt.plot(xax, eta, 'b-')
plt.xlabel("J, A/m^2")
plt.ylabel("eta")
plt.show()