import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

h = 6.626e-34
hbar = 1.051e-34
hbarel = 6.582e-16
e = 1.6e-19
V0 = 55e-3
da = 6e-9
mel = 9.11e-31
me = 0.067*mel
ml = 0.04*mel
mh = 0.3*mel
# def find_k(a,b, eps):
#     delt = np.inf
#     while delt < eps:

def func_SCH(k, args):
    m = args
    k0 = np.sqrt(2*m*e*V0)*1e-9/h
    if np.tan(k*da/2) <0:
        return abs(abs(np.sin(k*da/2))-k/k0)
    else:
        return abs(abs(np.cos(k * da / 2)) - k / k0)

res = minimize(func_SCH, np.array(0), tol = 1e-6, args = mel)
kel = res.x[0]
res = minimize(func_SCH, np.array(0), tol = 1e-6, args = ml)
kl = res.x[0]
res = minimize(func_SCH, np.array(0), tol = 1e-6, args = mh)
kh = res.x[0]

xax = np.linspace(0, 0.25,100)
plt.plot(xax, [func_SCH(x, mel) for x in xax], 'b-', label = 'електрон')
plt.plot(xax, [func_SCH(x, ml) for x in xax], 'r-', label = 'легка дірка')
plt.plot(xax, [func_SCH(x, mh) for x in xax], 'g-', label = 'важка дірка')
plt.xlabel('k, nm')
plt.ylabel('f')
plt.legend()
plt.show()

A2_f = lambda k: 1/np.sqrt(da/2 +1/k)
E_f = lambda k,m: hbar**2 * k**2 / (2 *m)/e
alph_f = lambda k,m: np.sqrt(2*m*e*E_f(k,m))/hbar

def psi(x, k, n,m):
    A2 = A2_f(k)
    alph = alph_f(k,m)
    #print(A2, alph)
    if n % 2 == 0:
        if x < - da/2:
            return -2 * A2 * np.exp(k*(da/2 +x)) * np.sin(alph *da/2)
        if abs(x)<= da/2:
            return 2 * A2 * np.sin(alph*x)
        if x > da/2:
            return 2 * A2 * np.exp(k * (da / 2 - x)) * np.sin(alph * da / 2)
    else:
        if x < - da/2:
            return 2 * A2 * np.exp(k*(da/2 +x)) * np.cos(alph *da/2)
        if abs(x) <= da/2:

            return 2 * A2 * np.cos(alph*x)
        if x > da/2:
            return 2 * A2 * np.exp(k * (da / 2 - x)) * np.cos(alph * da / 2)


xax = np.linspace(-2*da, 2*da, 1000)

plt.figure(2)
plt.vlines(-da/2,0, V0, color = 'k')
plt.vlines(da/2,0, V0, color = 'k')
plt.hlines(0 , -da/2, da/2, color = 'k')
plt.plot(xax, np.array([psi(x, kel*1e9, 0,mel) for x in xax])*1e-6 +0.03, 'b-')
plt.plot(xax, np.array([psi(x, kel*1e9, 1,mel) for x in xax])*1e-6 +0.03, 'r-')
#plt.plot(xax, [psi(x, kel*1e9, 2,mel) for x in xax], 'y-')
plt.show()

plt.figure(3)
plt.vlines(-da/2,0, V0, color = 'k')
plt.vlines(da/2,0, V0, color = 'k')
plt.hlines(0 , -da/2, da/2, color = 'k')
plt.plot(xax, np.array([psi(x, kh*1e9, 0,mh) for x in xax])*1e-6 +0.03, 'b-')
plt.plot(xax, np.array([psi(x, kh*1e9, 1,mh) for x in xax])*1e-6 +0.03, 'r-')
#plt.plot(xax, [psi(x, kel*1e9, 2,mel) for x in xax], 'y-')
plt.show()

plt.figure(4)
plt.vlines(-da/2,0, V0, color = 'k')
plt.vlines(da/2,0, V0, color = 'k')
plt.hlines(0 , -da/2, da/2, color = 'k')
plt.plot(xax, np.array([psi(x, kl*1e9, 0,ml) for x in xax])*1e-6 +0.03, 'b-')
plt.plot(xax, np.array([psi(x, kl*1e9, 1,ml) for x in xax])*1e-6 +0.03, 'r-')
#plt.plot(xax, [psi(x, kel*1e9, 2,mel) for x in xax], 'y-')
plt.show()