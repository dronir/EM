

t_vscatter = """
&params
mediumFileName    = "${mFilename}"
hsFilename        = "${oFilename}"
sDistribution     = "constant"
resThetaE         = $rTheta
gridResHorizontal = 400
gridResVertical   = 10
rho               = 0.5
rhoAllowedError   = 100.0
sampleSeed        = 0

nThetaIn	  = $nTheta
thetaIn           = "$theta"

nOrders           = $nOrders
nSamplesPerOrder  = "$nSamples 1"

rf_applyFields    = ${applyFields}
rf_nFieldsPerMed  = $nFields
rf_spectrumType   = "${sType}"
rf_P              = $H
rf_std            = $std

brdfType          = "${brdfType}"

nThreads          = $nThreads
/
"""

t_dscatter = """
&params
mediumFileName    = "${mFilename}"
outFilename       = "${oFilename}"
sDistribution     = "constant"
gridResHorizontal = 400
gridResVertical   = 10
rho               = 0.5
rhoAllowedError   = 100.0
sampleSeed        = 0

nOrders           = 1
nSamplesPerOrder  = "$nSamples 1"

rf_applyFields    = ${applyFields}
rf_nFieldsPerMed  = $nFields
rf_spectrumType   = "${sType}"
rf_P              = $H
rf_std            = $std

nThreads          = $nThreads
/
"""

t_submit = """
#!/bin/csh

#BSUB -L /bin/csh

#BSUB -J ${jobName}
#BSUB -e error.txt
#BSUB -o out.txt

#BSUB -M 300000
#BSUB -W 20:00
#BSUB -n 1

#BSUB -u "hannu.p.parviainen@helsinki.fi"
#BSUB -N
#BSUB -B

## run 
srun $$WRKDIR/vScatter/vScatter $$WRKDIR/vScatter/${brdfType}/$parFile

#bjobs -l $$LSB_JOBID
"""

