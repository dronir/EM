#!/usr/bin/python

import pylab as pl
import numpy as np

from numpy import abs, floor, sqrt, exp, log, cos, arccos, sin, arcsin, pi
from numpy.random import rand

from math import radians as rad, degrees as deg

from gather import *

HALF_PI     = 0.5*pi
TWO_PI      = 2.0*pi
FOUR_PI     = 4.0*pi
INV_FOUR_PI = 1.0 / (4.0 * pi)

rEps = 1e-7
iEps = 1e-24

def randExp():
    return -log(rand())

def H(w,x):
    g  = sqrt(1. - w)
    r0 = (1. - g) / (1. + g)
    return 1.0/(1.0 - w*x*( r0 + .5*(1. - 2.*r0*x) * log((1. + x) / x) ) )

def randCartHemiSphUniform(d):

    d[2] = rand()

    phi  = rand() * 2.0 * pi
    r    = sqrt(1.0 - d[2]**2)

    d[0] = r * cos(phi)
    d[1] = r * sin(phi)


def randCartSphUniform(d):

    d[2] = rand()

    phi  = rand() * 2.0 * pi
    r    = sqrt(1.0 - d[2]**2)

    d[0] = r * cos(phi)
    d[1] = r * sin(phi)

    if(rand() < 0.5):
        d[2] = -d[2]

def mcRadTransSS(i, w, mu, rTheta, nRays):
   
    hs = gth_hemisphere(rTheta,1,1)

    dTht   = HALF_PI / float(rTheta)
    dThtI  = float(rTheta) / HALF_PI

    nSplit    = 1
    nSplitInv = 1. / float(nSplit)

    for iRay in range(nRays):

        rI = cos(i)/float(nRays)
        rP = np.zeros(3)
        rD = np.array([sin(i), 0., -cos(i)])

        tau = randExp()

        ## Travel
        ##
        rP += rD * tau/mu

        ## Scattering
        ##
        for iSplit in range(nSplit):
            randCartHemiSphUniform(rD)

            rIt = nSplitInv * 0.5 * w * rI

            ## Attenuation
            ##
            l    = abs(rP[2] / rD[2])
            rIt *= exp(-l*mu)

            hs.addDataCar(rD, rIt)

    return hs 


def mcRadTrans(i, w, mu, rTheta, nRays):

    hs = gth_hemisphere(rTheta,1,1)
    
    dTht   = HALF_PI / float(rTheta)
    dThtI  = float(rTheta) / HALF_PI

    nSplit    = 1
    nSplitInv = 1. / float(nSplit)

    directWeight = 0.05

    rIo = cos(i) / float(nRays)
    rDo = np.array([sin(i), 0., -cos(i)])
    rPo = np.zeros(3)-0.01

    rDp = np.zeros(3)

    for iRay in range(nRays):

        rI = rIo
        rP = rPo.copy() 
        rD = rDo.copy() 

        while (rP[2] < 0.0) and (rI > iEps):
            tau = randExp()

            ## Travel
            ##
            rP += rD * tau/mu

            if(rP[2] < 0.0):

                ## Peeling
                ## 
                for iSplit in range(nSplit):
                    randCartHemiSphUniform(rDp)
                    l    = abs(rP[2] / rDp[2])       
         
                    hs.addDataCar(rDp, directWeight * nSplitInv * 0.5 * w * exp(-l*mu) * rI)

                ## Scattering
                ##
                randCartSphUniform(rD)
                rI *= (1. - directWeight) * w
            
        if(rI > iEps):
            hs.addDataCar(rD, rI)

    return hs
            
def mcSiRadTrans(i, e, w, mu, nRays):
    
    nSplit    = 1
    nSplitInv = 1. / float(nSplit)

    directWeight = 0.0

    rIo = 1.0 / float(nRays)
    rDo = np.array([sin(i), 0., -cos(i)])
    rPo = np.zeros(3)-0.01

    rDp = np.zeros(3)

    I = 0.0

    for iRay in range(nRays):

        rI = rIo
        rP = rPo.copy() 
        rD = rDo.copy() 

        while (rP[2] < 0.0) and (rI > iEps):
            tau = randExp()

            ## Travel
            ##
            rP += rD * tau/mu

            if(rP[2] < 0.0):

                rI *= w

                ## Peel
                ## 
                for iSplit in range(nSplit):
                    l   = abs(rP[2] / cos(e))  
                    I  += INV_FOUR_PI * exp(-l*mu) * rI
                   
                ## Scatter
                ##
                randCartSphUniform(rD)
                
            
    return I


def LS(i, e, w):
    return w * INV_FOUR_PI * cos(i) / (cos(i) + cos(e))

def CH(i, e, w):
    return w * INV_FOUR_PI * cos(i) / (cos(i) + cos(e)) * H(w, cos(i)) * H(w, cos(e))

## ANGLE INITIALIZATION
##
theta_i     = rad(65.0)
theta_e_res = 20
theta_e     = (np.arange(0,theta_e_res) + 0.5) * HALF_PI / theta_e_res

## MODEL INITIALIZATION
##
nRays = 15000
w     = 0.5
mu    = 1.0

## MODEL EVALUATION
##

MC   = mcRadTrans(theta_i, w, mu, theta_e_res, nRays)
MC.divideBySolidAngle()
 
I_mc = MC.toArray()
I_ls = LS(theta_i, theta_e, w)
I_ch = CH(theta_i, theta_e, w)

#te = rad(13.5)
#I_tt = cos(theta_i) * mcSiRadTrans(theta_i, te, w, mu, 3*nRays)/cos(te)

#print I_tt

#pl.plot([te],[I_tt], '+')

pl.plot(theta_e, I_mc.mean(1)/cos(theta_e))
pl.plot(theta_e, I_ls)
pl.plot(theta_e, I_ch,'o')
pl.show()
