#!/usr/bin/python2.4

import sys
from PyQt4 import Qt
from PyQt4 import QtGui
from PyQt4 import QtCore

import numpy      as np
import math

class xrGui(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)

        self.setWindowTitle('xrGui')


        ## ELEMENT LIST
        ##
        self.elementGrpB = QtGui.QGroupBox("Elements")
        self.elementGrpB.setMinimumHeight(250)
        self.elementGrpB.setMaximumWidth(200)

        self.elementList = QtGui.QListWidget(self)
        #self.elementList.addItems(self.hs.Elements)
        self.elementList.setCurrentRow(0)
        self.elementList.setEnabled(False)
        self.elementList.setMinimumHeight(250)
        self.elementList.setMaximumWidth(200)

        elementLayout = QtGui.QVBoxLayout()
        elementLayout.addWidget(self.elementList)
        
        self.elementGrpB.setLayout(elementLayout)

        ## THETA LIST
        ##
        self.thetaGrpB = QtGui.QGroupBox("Angle of incidence")
        self.thetaGrpB.setMinimumHeight(300)
        self.thetaGrpB.setMaximumWidth(200)

        self.thetaList = QtGui.QListWidget(self)
        #for i in range(self.hs.nThetaI):
        #    self.thetaList.addItem("%6.2f" % math.degrees(self.hs.thetaI[i]))
        self.thetaList.setCurrentRow(0)
        self.thetaList.setMinimumHeight(100)
        self.thetaList.setMaximumWidth(200)


        thetaLayout = QtGui.QVBoxLayout()
        thetaLayout.addWidget(self.thetaList)
        
        self.thetaGrpB.setLayout(thetaLayout)


        ## PLOT TYPE
        ##
        self.plotStyleGroup    = QtGui.QButtonGroup(self)
        self.plotStyleGroupBox = QtGui.QGroupBox("Plot type")
        self.plotStyleGroupBox.setMaximumWidth(200)

        pAbsolute = QtGui.QRadioButton("Absolute")
        pRelative = QtGui.QRadioButton("Relative")
        pNormalRe = QtGui.QRadioButton("Normalized relative")

        pAbsolute.setChecked(True)

        self.plotStyleGroup.addButton(pAbsolute, 0)
        self.plotStyleGroup.addButton(pRelative, 1)
        self.plotStyleGroup.addButton(pNormalRe, 2)

        plotTypeLayout = QtGui.QVBoxLayout()

        plotTypeLayout.addWidget(pAbsolute)
        plotTypeLayout.addWidget(pRelative)
        plotTypeLayout.addWidget(pNormalRe)

        self.plotStyleGroupBox.setLayout(plotTypeLayout)

        ## 
        ##

        
        ## QUIT BUTTON
        ##
        self.qb = QtGui.QPushButton("Quit")
        self.qb.setMaximumWidth(100)

        self.hb = QtGui.QPushButton("Toggle Plottype")
        self.hb.setMaximumWidth(150)

        self.toggleAn = QtGui.QCheckBox("Solution")
        self.toggleAn.setChecked(True)

        self.showInfo = QtGui.QCheckBox("Show info")
        self.showInfo.setChecked(True)

        self.showPlots = QtGui.QCheckBox("Show plots")
        self.showPlots.setChecked(True)

        ## PLOTS
        ## 
        #self.linePlots = linePlot(self.hs, self.data, self)
        #self.mu        = muPlot(self.hs, self)
        #self.spectrum  = spectrumPlot(self.hs, self)
        #self.hsPlots   = hsPlot(self.hs, self)

        #self.hsPlots.setVisible(False)

        #self.dPlot = Qwt.QwtPlot()
        #self.dPlot.setAxisScale(Qwt.QwtPlot.xBottom, data.energy.min(), data.energy.max())
        #self.dPlot.setAxisScale(Qwt.QwtPlot.yLeft, 0, data.muExt.max())
        #self.dPlot.setCanvasBackground(Qt.Qt.white)

        #dPlotData  = Qwt.QwtPlotCurve('Density Structure')
        #dPlotData.setPen(Qt.QPen(Qt.Qt.black))

        #dPlotData.setData(self.hs.medDensitymap[0,:,0], self.hs.medDensitymap[1,:,0])
        #dPlotData.attach(self.dPlot)


        ## SET LAYOUT
        ##
        sbox = QtGui.QHBoxLayout()

        rbox = QtGui.QVBoxLayout()

        hbox = QtGui.QVBoxLayout()
        hbox.addWidget(self.plotStyleGroupBox)
        hbox.addWidget(self.elementGrpB)
        hbox.addWidget(self.thetaGrpB)
        hbox.addStretch(1)
        hbox.addWidget(self.toggleAn)
        hbox.addWidget(self.showPlots)
        hbox.addWidget(self.showInfo)
        hbox.addWidget(self.hb)
        hbox.addWidget(self.qb)

        sbox.addLayout(hbox)
        sbox.addSpacing(50)
        sbox.addLayout(rbox)

        #sbox.addWidget(self.linePlots)
        #sbox.addWidget(self.hsPlots)

        
        infoLayout = QtGui.QVBoxLayout()
        #infoLayout.addWidget(self.mu)
        #if(self.hs.spcType=='Spectrum'):
        #    infoLayout.addWidget(self.spectrum)
        #infoLayout.addWidget(self.dPlot)

        sbox.addLayout(infoLayout)

        self.setLayout(sbox)

        self.resize(800, 1000)

        ## SET CONNECTIONS
        ##

        #self.connect(self.plotStyleGroup, QtCore.SIGNAL('buttonClicked(int)'), self.plot)

        #self.connect(self.elementList, QtCore.SIGNAL('itemSelectionChanged()'), self.plot)
        #self.connect(self.thetaList, QtCore.SIGNAL('itemSelectionChanged()'), self.plot)

        #self.connect(self.showPlots, QtCore.SIGNAL('stateChanged(int)'), self.togglePlots)
        #self.connect(self.showInfo, QtCore.SIGNAL('clicked()'), self.toggleInfo)
        #self.connect(self.qb, QtCore.SIGNAL('clicked()'), QtGui.qApp, QtCore.SLOT('quit()'))

        #self.connect(self.hb, QtCore.SIGNAL('clicked()'), self.hide)

    
xrGuiApp = QtGui.QApplication(sys.argv)
xrGuiWin = xrGui()

xrGuiWin.show()
xrGuiApp.exec_()
