import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

h = 6.626e-34
hbar = 1.051e-34
hbarel = 6.582e-16
e = 1.6e-19
V0 = 55e-3
da = 6#nm
me = 9.11e-31
mel = 0.067*me
ml = 0.04*me
mh = 0.3*me
# def find_k(a,b, eps):
#     delt = np.inf
#     while delt < eps:

def func_SCH(k, *args):
    m = args[0]
    k0 = np.sqrt(2*m*e*V0)*1e-9/hbar
    if np.tan(k*da/2) <0:
        return abs(abs(np.sin(k*da/2))-k/k0)
    else:
        return abs(abs(np.cos(k * da / 2)) - k / k0)

res = minimize(func_SCH, np.array(0), tol = 1e-6, args = tuple([mel]))
kel = res.x[0]
res = minimize(func_SCH, np.array(0), tol = 1e-6, args = tuple([ml]))
kl = res.x[0]
res = minimize(func_SCH, np.array(0), tol = 1e-6, args = tuple([mh]))
kh = res.x[0]

res = minimize(func_SCH, np.array(0.6), tol = 1e-6, args = tuple([mel]),bounds = ((0.5,1),))
kel2 = res.x[0]
res = minimize(func_SCH, np.array(0.6), tol = 1e-6, args = tuple([ml]),bounds = ((0.5,1),) )
kl2 = res.x[0]
res = minimize(func_SCH, np.array(0.6), tol = 1e-6, args = tuple([mh]),bounds = ((0.6,1),))
kh2 = res.x[0]

res = minimize(func_SCH, np.array(1.1), tol = 1e-6, args = tuple([mel]),bounds = ((1,2),))
kel3 = res.x[0]
res = minimize(func_SCH, np.array(1.1), tol = 1e-6, args = tuple([ml]),bounds = ((1,2),) )
kl3 = res.x[0]
res = minimize(func_SCH, np.array(1.1), tol = 1e-6, args = tuple([mh]),bounds = ((1,2),))
kh3 = res.x[0]

xax = np.linspace(0, 1.5,100)
plt.plot(xax, [func_SCH(x, mel) for x in xax], 'b-', label = 'електрон')
plt.plot(xax, [func_SCH(x, ml) for x in xax], 'r-', label = 'легка дірка')
plt.plot(xax, [func_SCH(x, mh) for x in xax], 'g-', label = 'важка дірка')
plt.xlabel('k, nm')
plt.ylabel('f')
plt.legend()
plt.show()

A2_f = lambda k: 1/np.sqrt(da/2 +1/k)
E_f = lambda k,m: hbar**2 * (k*1e9)**2 / (2 *m)/e
alph_f = lambda k,m: np.sqrt(2*m*e*E_f(k,m))/hbar*1e-9

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


# plt.figure(2)
# plt.vlines(-da/2,0, V0, color = 'k')
# plt.vlines(da/2,0, V0, color = 'k')
# plt.hlines(0 , -da/2, da/2, color = 'k')
# plt.plot(xax, np.array([psi(x, kel, 0,mel)*1e-3 + 0.03 for x in xax]), 'b-')
# plt.plot(xax, np.array([psi(x, kel, 1,mel)*1e-3 + 0.03 for x in xax]), 'r-')
# plt.plot(xax, np.array([psi(x, kel2, 0,mel)*1e-3 + 0.02 for x in xax]), 'b--')
# plt.plot(xax, np.array([psi(x, kel2, 1,mel)*1e-3 + 0.02 for x in xax]), 'r--')
# plt.plot(xax, np.array([psi(x, kel3, 0,mel)*1e-3 + 0.01 for x in xax]), 'b.')
# plt.plot(xax, np.array([psi(x, kel3, 1,mel)*1e-3 + 0.01 for x in xax]), 'r.')
# #plt.plot(xax, [psi(x, kel*1e9, 2,mel) for x in xax], 'y-')
# plt.show()
#
# plt.figure(3)
# plt.vlines(-da/2,0, V0, color = 'k')
# plt.vlines(da/2,0, V0, color = 'k')
# plt.hlines(0 , -da/2, da/2, color = 'k')
# plt.plot(xax, np.array([psi(x, kh, 0,mh)*1e-3 + 0.03 for x in xax]), 'b-')
# plt.plot(xax, np.array([psi(x, kh, 1,mh)*1e-3 + 0.03 for x in xax]), 'r-')
# plt.plot(xax, np.array([psi(x, kh2, 0,mh)*1e-3 + 0.02 for x in xax]), 'b--')
# plt.plot(xax, np.array([psi(x, kh2, 1,mh)*1e-3 + 0.02 for x in xax]), 'r--')
# plt.plot(xax, np.array([psi(x, kh3, 0,mh)*1e-3 + 0.01 for x in xax]), 'b.')
# plt.plot(xax, np.array([psi(x, kh3, 1,mh)*1e-3 + 0.01 for x in xax]), 'r.')
# plt.xlabel("x, nm")
# #plt.plot(xax, [psi(x, kel*1e9, 2,mel) for x in xax], 'y-')
# plt.show()
#
# plt.figure(4)
# plt.vlines(-da/2,0, V0, color = 'k')
# plt.vlines(da/2,0, V0, color = 'k')
# plt.hlines(0 , -da/2, da/2, color = 'k')
# plt.plot(xax, np.array([psi(x, kl, 0,ml)*1e-3 + 0.03 for x in xax]), 'b-')
# plt.plot(xax, np.array([psi(x, kl, 1,ml)*1e-3 + 0.03 for x in xax]), 'r-')
# plt.plot(xax, np.array([psi(x, kl2, 0,ml)*1e-3 + 0.02 for x in xax]), 'b--')
# plt.plot(xax, np.array([psi(x, kl2, 1,ml)*1e-3 + 0.02 for x in xax]), 'r--')
# plt.plot(xax, np.array([psi(x, kl3, 0,ml)*1e-3 + 0.01 for x in xax]), 'b.')
# plt.plot(xax, np.array([psi(x, kl3, 1,ml)*1e-3 + 0.01 for x in xax]), 'r.')
# #plt.plot(xax, [psi(x, kel*1e9, 2,mel) for x in xax], 'y-')
# plt.show()

def draw(karr, m, fig_num):
    plt.figure(fig_num)
    plt.vlines(-da/2,0, V0, color = 'k')
    plt.vlines(da/2,0, V0, color = 'k')
    plt.hlines(0 , -da/2, da/2, color = 'k')
    pointer = ['b-','r-', 'b--','r--','b.','r.']
    for i in range(len(karr)):
        plt.plot(xax, np.array([psi(x, karr[i], 0,m)*2e-3 + 0.04 - i*1e-2 for x in xax]), pointer[2*i])
        plt.plot(xax, np.array([psi(x, karr[i], 1,m)*2e-3 + 0.04 - i*1e-2 for x in xax]), pointer[2*i +1])
    #plt.plot(xax, [psi(x, kel*1e9, 2,mel) for x in xax], 'y-')
    plt.xlabel("x, nm")
    plt.show()

draw([kel,kel2,kel3], mel, 2)
draw([kl,kl2,kl3], ml, 3)
draw([kh,kh2,kh3], mh, 4)


karr = [kel, kl,kh,kel2,kl2,kh2,kel3,kl3,kh3]
marr = [mel,ml,mh,mel,ml,mh,mel,ml,mh]
Earr = [E_f(karr[i],marr[i]) for i in range(9)]
A2arr = [A2_f(karr[i]) for i in range(9)]
print(np.array(np.transpose([Earr,A2arr,karr])))