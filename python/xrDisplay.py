#!/usr/bin/python

import sys
from PyQt4 import Qt
from PyQt4 import QtGui
from PyQt4 import QtCore

import PyQt4.Qwt5 as Qwt
import numpy      as np

import math

from gather import *

class hsData(Qwt.QwtRasterData):

    def __init__(self, data, min, max):
        Qwt.QwtRasterData.__init__(self, Qt.QRectF(0.0, 0.0, 360.0, 90.0))

        self.data   = np.flipud(data)

        self.dTheta = 1.0 / (90.001  / (data.shape[0]))
        self.dPhi   = 1.0 / (360.001 / (data.shape[1]))

        self.min = min
        self.max = max

    def copy(self):
        return self

    def range(self):
        #return Qwt.QwtDoubleInterval(self.min, self.max);
        return Qwt.QwtDoubleInterval(self.data.min(), self.data.max());

    def value(self, phi, theta):
        return self.data[self.dTheta*theta, self.dPhi*phi]


class heightmapData(Qwt.QwtRasterData):
    def __init__(self, data, i):
        Qwt.QwtRasterData.__init__(self, Qt.QRectF(0.0, 0.0, data.shape[0]-1, data.shape[1]-1))

        self.i    = i
        self.data = data
        self.dx   = 1.0

    def copy(self):
        return self

    def range(self):
        return Qwt.QwtDoubleInterval(self.data.min(), self.data.max());


    def value(self, x, y):
        return self.data[self.i, x, y]



class hsPlot(QtGui.QWidget):
    def __init__(self, hs, iTheta=0, parent=None):
        QtGui.QWidget.__init__(self, parent)

        layOut = QtGui.QVBoxLayout()

        self.hs    = hs

        self.plots = []
        self.specs = []

        colorMap = Qwt.QwtLinearColorMap(Qt.Qt.black, Qt.Qt.white)

        for i in range(len(hs.Elements)):
            self.plots.append(Qwt.QwtPlot())
            self.specs.append(Qwt.QwtPlotSpectrogram())

            self.plots[i].setCanvasBackground(Qt.Qt.white)

            self.plots[i].setAxisScale(Qwt.QwtPlot.xBottom, 0.0, 360.0)
            self.plots[i].enableAxis(Qwt.QwtPlot.xBottom, False)

            self.plots[i].setAxisScale(Qwt.QwtPlot.yLeft, 0.0, 90.0)
            self.plots[i].enableAxis(Qwt.QwtPlot.yLeft, False)


        for i in range(len(hs.Elements)):
            data = hs.toArray(0, i*2) + hs.toArray(0, i*2+1)
            
            self.specs[i].setData(hsData(data, hs.data.min(), hs.data.max()))
            self.specs[i].setColorMap(colorMap)
            self.specs[i].attach(self.plots[i])

            layOut.addWidget(self.plots[i])

        self.setBaseSize(400,100*len(hs.Elements))
        self.setMinimumSize(200,50*len(hs.Elements))

        self.setLayout(layOut)

    def replot(self, pType, iRefElem, iTht):
        
        for i in range(len(self.hs.Elements)):

            if(pType == 0):
                data = self.hs.toArray(iTht, i*2) + self.hs.toArray(iTht, i*2+1)

            elif(pType == 1):
                data = (self.hs.toArray(iTht, i*2) + self.hs.toArray(iTht, i*2+1)) / (self.hs.toArray(iTht, iRefElem*2) + self.hs.toArray(iTht, iRefElem*2+1))

            elif(pType == 2):
                data = (self.hs.toArray(iTht, i*2) + self.hs.toArray(iTht, i*2+1)) / (self.hs.toArray(iTht, iRefElem*2) + self.hs.toArray(iTht, iRefElem*2+1))

                data /= data[0,0]

            #data = self.hs.toArray(iTht, i*2) + self.hs.toArray(iTht, i*2+1)
            
            self.specs[i].setData(hsData(data,  self.hs.data.min(), self.hs.data.max()))
            self.plots[i].replot()


class heightmapPlot(QtGui.QWidget):
    def __init__(self, hs, iTheta=0, parent=None):
        QtGui.QWidget.__init__(self, parent)

        layOut = QtGui.QVBoxLayout()

        self.hs    = hs

        self.plots = []
        self.specs = []

        colorMap = Qwt.QwtLinearColorMap(Qt.Qt.black, Qt.Qt.white)

        for i in range(hs.medHeightmap.shape[2]):
            self.plots.append(Qwt.QwtPlot())
            self.specs.append(Qwt.QwtPlotSpectrogram())

            self.plots[i].setCanvasBackground(Qt.Qt.white)

            #self.plots[i].setAxisScale(Qwt.QwtPlot.xBottom, 0.0, 360.0)
            self.plots[i].enableAxis(Qwt.QwtPlot.xBottom, False)

            #self.plots[i].setAxisScale(Qwt.QwtPlot.yLeft, 0.0, 90.0)
            self.plots[i].enableAxis(Qwt.QwtPlot.yLeft, False)


        for i in range(hs.medHeightmap.shape[2]):
            
            self.specs[i].setData(heightmapData(hs.medHeightmap[:,:,i]))
            self.specs[i].setColorMap(colorMap)
            self.specs[i].attach(self.plots[i])

            layOut.addWidget(self.plots[i])

        #self.setBaseSize(400,100*len(hs.Elements))
        self.setMinimumSize(200,200*hs.medHeightmap.shape[2])

        self.setLayout(layOut)


if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)

    hs = xrHemisphere() 
    hs.load(sys.argv[1])

    if(hs.Method == "MonteCarlo"):
        hs.divideBySolidAngle()

    #icon = hsPlot(hs)
    icon = heightmapPlot(hs)

    icon.show()
    app.exec_()
