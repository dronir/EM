#!/usr/bin/python

import sys
#from PyQt4 import Qt
#from PyQt4 import QtGui
#from PyQt4 import QtCore

import qt
import Qwt5 as Qwt

import numpy      as np

import math

from gather import *

class elementList(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)

        self.setWindowTitle('Element List')

        self.pt = QtGui.QCheckBox("Relative plots", self)

        self.nt = QtGui.QCheckBox("Normalize plots", self)
        self.nt.setEnabled(False)

        self.el = QtGui.QListWidget(self)
        self.el.addItems(hs.Elements)
        self.el.setCurrentRow(0)
        self.el.setEnabled(False)
        self.el.setMinimumHeight(350)
        self.el.setMaximumWidth(150)

        self.tl = QtGui.QListWidget(self)
        for i in range(hs.nThetaI):
            self.tl.addItem("%6.2f" % math.degrees(hs.thetaI[i]))
        self.tl.setCurrentRow(0)
        self.tl.setMinimumHeight(150)
        self.tl.setMaximumWidth(150)


        self.qb = QtGui.QPushButton("Quit")
        self.qb.setMaximumWidth(100)

        ## SET PLOT ATTRIBUTES
        ## 
        self.ratioPlots = []
        for i in range(len(hs.Elements)):
            self.ratioPlots.append( Qwt.QwtPlot(self))
            self.ratioPlots[i].enableAxis(Qwt.QwtPlot.xBottom, False)

            self.ratioPlots[i].setCanvasBackground(Qt.Qt.white)

            #self.ratioPlots[i].plotLayout().setCanvasMargin(0)
            #self.ratioPlots[i].plotLayout().setAlignCanvasToScales(True)

            self.ratioPlots[i].setAxisScale(Qwt.QwtPlot.xBottom, 0, 90)
            self.ratioPlots[i].setAxisMaxMajor(Qwt.QwtPlot.yLeft, 4)
            #self.ratioPlots[i].axisWidget(Qwt.QwtPlot.yLeft).setBorderDist(50,60)

        ## LOAD DATA
        ##
        self.data     = []
        for iTht in range(hs.nThetaI):
            self.data.append([])
            for iElem in range(len(hs.Elements)):
                self.data[iTht].append((hs.toArray(set=iTht, lvl=iElem*2) + hs.toArray(set=iTht, lvl=iElem*2+1)).mean(1))


        ## PLOT
        ##
        self.plotData = []
        x = np.linspace(0, 90, hs.resTheta)
        for iElem in range(len(hs.Elements)):
            self.plotData.append(Qwt.QwtPlotCurve('y = sin(x)'))
            self.plotData[iElem].setPen(Qt.QPen(Qt.Qt.red))

            y = self.data[0][iElem]

            self.plotData[iElem].setData(x, y)
            self.plotData[iElem].attach(self.ratioPlots[iElem])


        ## SET LAYOUT
        ##
        sbox = QtGui.QHBoxLayout()
        rbox = QtGui.QVBoxLayout()

        hbox = QtGui.QVBoxLayout()
        hbox.addWidget(self.el)
        hbox.addWidget(self.pt)
        hbox.addWidget(self.nt)
        hbox.addSpacing(50)
        hbox.addWidget(self.tl)
        hbox.addStretch(1)
        hbox.addWidget(self.qb)

        for i in range(len(hs.Elements)):
            rbox.addWidget(self.ratioPlots[i])

        sbox.addLayout(hbox)
        sbox.addSpacing(50)
        sbox.addLayout(rbox)

        self.setLayout(sbox)

        self.resize(800, 1000)

        ## SET CONNECTIONS
        ##
        self.connect(self.el, QtCore.SIGNAL('itemSelectionChanged()'), self.plot)
        self.connect(self.tl, QtCore.SIGNAL('itemSelectionChanged()'), self.plot)

        self.connect(self.pt, QtCore.SIGNAL('stateChanged(int)'), self.changeRel)
        self.connect(self.nt, QtCore.SIGNAL('stateChanged(int)'), self.changeRel)

        self.connect(self.qb, QtCore.SIGNAL('clicked()'), QtGui.qApp, QtCore.SLOT('quit()'))


    def plot(self):

        iTht = self.tl.currentRow()

        x = np.linspace(0, 90, hs.resTheta)
        for iElem in range(len(hs.Elements)):

            if(self.pt.isChecked()):
                y = self.data[iTht][iElem] / self.data[iTht][self.el.currentRow()]

                if(self.nt.isChecked()):
                    y /= y[0]

            else:
                y = self.data[iTht][iElem]

            self.plotData[iElem].setData(x, y)
            self.plotData[iElem].attach(self.ratioPlots[iElem])

            self.ratioPlots[iElem].replot()

    def changeRel(self):
        
        self.nt.setEnabled(self.pt.isChecked()) 
        self.el.setEnabled(self.pt.isChecked())
        self.plot()


app = QtGui.QApplication(sys.argv)

hs = xrHemisphere() 
hs.load(sys.argv[1])
hs.divideBySolidAngle()

icon = elementList()

icon.show()
app.exec_()
