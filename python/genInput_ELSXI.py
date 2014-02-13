#!/usr/bin/python

import sys
import os

from string import Template

medRoot = sys.argv[1]

rt = 180
ns = 1000
nf = 10
nOrders = 1

nt = 9

t = []
for i in range(nt):
    t.append(float(i)*10.0)

nThreads = 1

subTemplateF = open('/home/hannu/work/Projects/BepiColombo/XR/python/template_submit.txt','r')
parTemplateF = open('/home/hannu/work/Projects/BepiColombo/XR/python/template_sc.txt','r')

subTemplate = Template(subTemplateF.read())
parTemplate = Template(parTemplateF.read())

brdfType    = sys.argv[2]

sTypes = ['smooth','Gaussian','fBm']
sTypes = ['smooth']
HVals  = [0.4, 0.6, 0.8]
lVals  = [0.5, 1.0, 1.5]

SVals  = [0.03, 0.04, 0.05, 0.06]

rVals  = ['020','035','050']

oVals  = [0.1, 0.5, 0.9]

for sType in sTypes:
    if sType == 'smooth':
        for rho in rVals:
            for omega in oVals:

                medFilename = medRoot+'/medium_rho_'+rho+'.nc'

                rhoStr = 'R'+rho 
                oStr   = 'w%04.2f' %omega

                subFilename = 'scatter_smooth_' + rhoStr + '_' + oStr + '.sub'
                parFilename = 'scatter_smooth_' + rhoStr + '_' + oStr + '.in'
                outFilename = brdfType+'_smooth_' + rhoStr  + '_' + oStr + '.vsc'

                f = open(parFilename, 'w')

                f.write(parTemplate.substitute(mFilename=medFilename,
                                               oFilename=outFilename,
                                               rTheta=rt,
                                               nOrders=nOrders,
                                               nSamples=ns*nf,
                                               nTheta=nt,
                                               theta=str(t).strip('[]'),
                                               nFields=1,
                                               H=0.0,
                                               std=0.0,
                                               omega = omega,
                                               mu=0.0,
                                               applyFields = '.false.',
                                               sType = 'fBm',
                                               nThreads = nThreads,
                                               brdfType = brdfType))

                f.close()

                f = open(subFilename, 'w')
                f.write(subTemplate.substitute(jobName='vScat'+rhoStr, parFile=parFilename, brdfType=brdfType))
                f.close()


    else:
        appFields = '.true.'
   
        if sType == 'Gaussian':
            P1 = lVals
        else:
            P1 = HVals

        for P in P1:
            for S in SVals:
                for rho in rVals:

                    medFilename = medRoot+'/medium_rho_'+rho+'.nc'

                    rhoStr = 'R'+rho 
                    hStr = '__P' + ('%4.2f' %P)#.replace('.','_')#[-1]
                    sStr = '__S' + ('%4.2f' %S)#[-2:]

                    subFilename = 'scatter_' + sType + '__' + rhoStr + hStr + sStr + '.sub'
                    parFilename = 'scatter_' + sType + '__' + rhoStr + hStr + sStr + '.in'
                    outFilename = 'LommelSeeliger_' + sType + '__' + rhoStr + hStr + sStr + '.vsc'
                        

                    f = open(parFilename, 'w')
                    
                    f.write(parTemplate.substitute(mFilename=medFilename,
                                                   oFilename=outFilename,
                                                   rTheta=rt,
                                                   nOrders=nOrders,
                                                   nSamples=ns,
                                                   nTheta=nt,
                                                   theta=str(t).strip('[]'),
                                                   nFields=nf,
                                                   H=P,
                                                   std=S,
                                                   applyFields = '.true.',
                                                   sType = sType,
                                                   nThreads = nThreads))
                    f.close()


                    f = open(subFilename, 'w')
                    f.write(subTemplate.substitute(jobName='vScat'+rhoStr+hStr+sStr,
                                                   parFile=parFilename))
                    f.close()
