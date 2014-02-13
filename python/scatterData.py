import Scientific.IO.NetCDF as nc
import numpy as np
import sys
import math
import pylab as pl
import matplotlib.colors as colors

class scatterData:
  """Class implementing the gathering hemisphere."""
  resTheta    = 0
  dTheta      = 0.0
  dThetaInv   = 0.0
  rows        = []
  resPhi      = []
  dPhi        = []
  dPhiInv     = []

  dSetName    = []

  rghI        = 0
  stdI        = 0

  def load(self, fName):
    """
    Loads the scattering data from a netCDF file.

    Returns: nothing
    """

    try:
      dFile = nc.NetCDFFile(fName, "r")
    except IOError:
      print "Error reading file, exiting."
      sys.exit()

    self.data   = np.array(dFile.variables['scatterData'].getValue()) * 8.0 * np.pi
    self.thetaE = np.array(dFile.variables['thetaE'].getValue())
    self.thetaI = np.array(dFile.variables['thetaI'].getValue())
    self.srfStd = np.array(dFile.variables['std'].getValue())
    self.roughP = np.array(dFile.variables['roughP'].getValue())

    self.rho    = float(dFile.variables['scatterData'].volumeDensity)

    self.dTe    = np.abs(self.thetaE[1] - self.thetaE[0])
    self.dTi    = np.abs(self.thetaI[1] - self.thetaI[0])
    self.dStd   = np.abs(self.srfStd[1] - self.srfStd[0])
    self.dRgh   = np.abs(self.roughP[1] - self.roughP[0])

    dFile.close()

  def evalCont(self, rP, std, tE, tI):
    """
    
    """

    rI = int((rP  - self.roughP[0]) * self.dRgh)
    sI = int((std - self.srfStd[0]) * self.dStd)
    iI = int((tI  - self.thetaI[0]) * self.dTi)
    eI = int((tE  - self.thetaE[0]) * self.dTe)

    return self.data[rI, sI, eI, iI]

  def evalDisc(self, rI, sI, tE, tI):
    """
    
    """

    iI = int((tI  - self.thetaI[0]) * self.dTi)
    eI = int((tE  - self.thetaE[0]) * self.dTe)

    return self.data[rI, sI, eI, iI]
