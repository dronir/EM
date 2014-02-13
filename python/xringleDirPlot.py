import Scientific.IO.NetCDF as nc
import numpy as np
import sys
import math
import pylab as pl
import matplotlib.colors as colors

from numpy import floor, sqrt, sin, cos, arccos, arctan2, pi


class gth_hemisphere:
  """Class implementing the gathering hemisphere."""

  def __init__(self, resTheta=1, nThetaI=1, nDataLevels=1):

    self.Type      = "Hemisphere"
    
    self.resTheta  = resTheta
    self.dTheta    = 0.5 * pi / float(resTheta)
    self.dThetaInv = 1.0 / self.dTheta

    self.nCells    = 1
    self.nThetaI   = nThetaI
    self.nLevels   = nDataLevels
    self.type      = type

    self.dPhi      = np.zeros(resTheta)
    self.dPhiInv   = np.zeros(resTheta)
    self.dA        = np.zeros(resTheta)
    self.mTheta    = np.zeros(resTheta)
    self.nPhi      = np.zeros(resTheta, np.int64)
    self.cIdx      = np.zeros(resTheta)
 
    self.phiRange  = 2.0 * pi

    dA0 = self.phiRange * (1.0 - cos(self.dTheta))

    self.nPhi    [0] = 1
    self.cIdx    [0] = 0
    self.dPhi    [0] = self.phiRange
    self.dPhiInv [0] = 1.0 / self.phiRange
    self.dA      [0] = dA0
    self.mTheta  [0] = 0.5 * self.dTheta

    for i in range(1, resTheta):
       dPhi = dA0 / (cos(i * self.dTheta) - cos((i+1) * self.dTheta))
       rPhi = round(self.phiRange / dPhi)
       dPhi = self.phiRange / float(rPhi)

       self.nPhi  [i] = rPhi
       self.dPhi    [i] = dPhi
       self.dPhiInv [i] = 1.0 / dPhi
       self.dA      [i] = dPhi * (cos(i * self.dTheta) - cos((i+1) * self.dTheta))
       self.mTheta  [i] = self.dTheta * (float(i) - 0.5)

       self.cIdx    [i] = self.cIdx[i-1] + self.nPhi[i-1]

       self.nCells = self.nCells + rPhi   

    self.dAMean = self.phiRange / float(self.nCells)

    self.data   = np.zeros([nThetaI, self.nCells, nDataLevels])
    self.weight = np.zeros([nThetaI, self.nCells, nDataLevels])

  def load(self, fName):
    """
    Loads the hemisphere data from a netCDF file.

    Returns: nothing
    """

    try:
      dFile = nc.NetCDFFile(fName, "r")
    except IOError:
      print "Error reading file, exiting."
      sys.exit()

    if "Hemisphere" not in dFile.variables.keys():
      print "Error: not a proper hemisphere file."
      sys.exit()


    if "Elements" in  dir(dFile):
      self.Elements = str(dFile.Elements).split()

    self.Type = dFile.Type

    self.nPhi      = np.array(dFile.nPhi)
    self.cIdx      = np.array(dFile.cIdx)

    ## Convert Fortran indices to numpy indices
    if self.cIdx[0] == 1:
      self.cIdx -= 1

    self.dPhi      = np.array(dFile.dPhi)
    self.dPhiInv   = 1.0 / self.dPhi  

    self.nThetaI   = int(dFile.nThetaI)
    self.nLevels   = int(dFile.nLevels)
    self.resTheta  = int(dFile.nThetaE)

    self.dTheta    = 0.5 * math.pi / float(self.resTheta)
    self.dThetaInv = 1.0/self.dTheta

    self.dA        = dFile.dA

    self.data      = np.array(dFile.variables['Hemisphere'].getValue())

    dFile.close()

  def divideBySolidAngle(self):
    
    for i in range(self.resTheta):
      self.data[:, self.cIdx[i] : self.cIdx[i] + self.nPhi[i], :] /= self.dA[i]

  def carDirToCell(self, D):

    r     = sqrt    ( (D**2).sum()       )
    theta = arccos  ( D[2] / r           )
    phi   = arctan2 ( D[1] / r, D[0] / r )
    
    if( phi < 0.0 ):
      phi = 2.0*pi + phi
    
    t = floor( theta * self.dThetaInv  )
    p = floor( phi   * self.dPhiInv[t] )
    
    return self.cIdx[t] + p


  def addDataCar(self, D, v, set=0, lvl=0):

    c = self.carDirToCell(D)
    self.data[set, c, lvl] += v

    
  def toArray(self, set=0, lvl=0):
    """
    Unpacks the gathering hemisphere into a 2-dimensional array.

    Returns: numpy.array
    """

    resTheta = self.resTheta
    resPhi   = self.nPhi.max()

    if(self.Type == 'Hemisphere'):
      dp       = math.pi * 2.0 / float(resPhi)
    else:
      dp       = math.pi / float(resPhi)
          
    data = np.zeros([resTheta, resPhi])   

    #print self.data[lvl,:,set]

    for i in range(resTheta):
      dPhiI = dp * self.dPhiInv[i]
      for j in range(resPhi):
         data[i,j] = self.data[lvl, self.cIdx[i] + int(math.floor(j * dPhiI)), set]

    return data


  def phiSlice(self, theta, set=0, lvl=0):
    """

    Returns: numpy.array
    """
 
    iTheta = int(math.floor(theta * self.dThetaInv))
    resPhi = self.nPhi[iTheta]
    dPhi   = self.dPhi[iTheta]

    data   = np.zeros([resPhi,2])

    for i in range(resPhi):
      data[i,0] = (i + 0.5) * dPhi
      data[i,1] = self.rows[set][iTheta][i,lvl]

    return data

  def thetaSlice(self, phi, set=0, lvl=0):
    """

    Returns: numpy.array
    """

    data = np.zeros([self.resTheta, 2])

    for i in range(self.resTheta):
      data[i,0] = (i+0.5) * self.dTheta
      data[i,1] = self.rows[set][i][phi * self.dPhiInv[i],lvl]

    return data

  def eval(self, thtI, thtE, phi):
    
    iThtI = int(math.floor(thtI))
    iThtE = int(math.floor(math.radians(thtE) * self.dThetaInv))
    iPhi  = int(math.floor(math.radians(phi) * self.dPhiInv[iThtE]))

    return self.rows[iThtI][iThtE][iPhi,0] * 4.0 * math.pi


class xrHemisphere(gth_hemisphere):
  
  def __init__(self, resTheta=1, nThetaI=1, nDataLevels=1):
    gth_hemisphere.__init__(self, resTheta, nThetaI, nDataLevels)


  def load(self, fName):
    """
    Loads the hemisphere data from a netCDF file.

    Returns: nothing
    """

    try:
      dFile = nc.NetCDFFile(fName, "r")
    except IOError:
      print "Error reading file, exiting."
      sys.exit()

    if "Hemisphere" not in dFile.variables.keys():
      print "Error: not a proper hemisphere file."
      sys.exit()

    try:
      self.Elements       = str(dFile.Elements).split()
      self.muAbs          = np.array(dFile.variables['muAbs'].getValue())
      self.muAbsCDF       = np.array(dFile.variables['muAbsCdf'].getValue())
      self.muExt          = np.array(dFile.variables['muExt'].getValue())
      self.Spectrum       = np.array(dFile.variables['Spectrum'].getValue())
      self.SpectrumCdf    = np.array(dFile.variables['SpectrumCdf'].getValue())
      self.SpectrumCdfInv = np.array(dFile.variables['SpectrumCdfInv'].getValue())

      self.SpectrumTest   = np.array(dFile.variables['SpectrumTest'].getValue())

    except (KeyError, AttributeError):
      print "Error: Malformed input file, missing data."
      sys.exit(1)

    self.Type = dFile.Type

    self.nPhi      = np.array(dFile.nPhi)
    self.cIdx      = np.array(dFile.cIdx)

    ## Convert Fortran indices to numpy indices
    if self.cIdx[0] == 1:
      self.cIdx -= 1

    self.dPhi      = np.array(dFile.dPhi)
    self.dPhiInv   = 1.0 / self.dPhi  

    self.nThetaI   = int(dFile.nThetaI)
    self.nLevels   = int(dFile.nLevels)
    self.resTheta  = int(dFile.nThetaE)

    self.dTheta    = 0.5 * math.pi / float(self.resTheta)
    self.dThetaInv = 1.0/self.dTheta

    self.dA        = dFile.dA

    self.data      = np.array(dFile.variables['Hemisphere'].getValue())

    dFile.close()

