#!/usr/bin/python
import pylab as pl


HS = (1.0 / (4.0 * pl.pi))

def computeMu(rho, A, E):
    r_e = 2.8179402894e-15
    m_u = 1.6605387820e-27
    h   = 4.1356673300e-15
    c   = 299792458.0
    
    n_a = (rho*1e3) / (A * m_u)
    lam = h * c / E

    return 2.0 * r_e * n_a * lam


def numModel(ti, te, mu1, mu2, n, dl):
    l = dl * pl.cos(ti)
    s = mu1 * dl
    
    u = pl.exp(-s)
    v = 1.0 - u

    a = mu2 * dl * pl.cos(ti) / pl.cos(te)

    I = 1.0
    S = pl.zeros(pl.size(te), pl.Float64)

    try:
        for i in range(1,n+1):
            A  = I * v
            S += A * pl.exp(-a * i)
            I *= u
    except:
        print "Halt dammit: ", i

    return HS * S


def LS(theta_i, theta_e):
    return HS / (pl.cos(theta_i) + pl.cos(theta_e))


def LS_HP1(theta_i, theta_e, mu1, mu2):
    return HS * mu1 / (mu2*pl.cos(theta_i) + mu1*pl.cos(theta_e))

def LS_HP2(theta_i, theta_e, mu1, mu2):
    return HS / ((mu2/mu1)*pl.cos(theta_i) + pl.cos(theta_e))


def LS_KM(theta_i, theta_e):
    return HS / (pl.cos(theta_i) + pl.cos(theta_e)) # * (1.0 + k * a)


print "Mu:", computeMu(7.874, 55.85, 7300.0) * 3.70552


a = pl.load("xr_inflsab.dat")

th_i = pl.pi * 0.250
th_e = a[:,0]

mu1   = 1.3
mu2   = 1.4

mu2   = 115121.238297125

muMin = 1.0
muMax = 1.35
muN   = 10

muArr = pl.load('FeMu.dat')
#muArr = pl.linspace(muMin, muMax, muN)

mu1   = muArr[5]

#print muArr

I  = pl.zeros(pl.size(th_e), pl.Float)

for mu in muArr:
    I += LS_HP2(th_i, th_e, mu, mu2) / pl.size(muArr)


print mu

I_XRR  = a[:,1]
I_ORI  = LS(th_i, th_e)
I_HP1  = LS_HP1(th_i, th_e, mu1, mu2)
I_HP2  = LS_HP2(th_i, th_e, mu1, mu2)
I_KM   = LS_KM(th_i, th_e)
I_INT  = numModel(th_i, th_e, mu1, mu2, 10000, 0.001)


pl.plot(th_e, I, c='blue')
pl.plot(th_e, I_KM * 1.0, c='black')
pl.plot(th_e, I_HP2, c='red')

#pl.plot(th_e, I_INT / pl.cos(th_e) )
#pl.plot(th_e, I_XRR / pl.cos(th_e) )

pl.show()
