#!/usr/bin/python

import sys
from PyQt4 import Qt
from PyQt4 import QtGui
from PyQt4 import QtCore

import PyQt4.Qwt5 as Qwt
import numpy      as np

import math

from xrDisplay import *

from gather import *

class spectrumPlot(QtGui.QWidget):
    def __init__(self, data, parent=None):
        QtGui.QWidget.__init__(self, parent)

        self.plot = Qwt.QwtPlot()
        self.plot.setCanvasBackground(Qt.Qt.white)

        self.plot.setAxisScale(Qwt.QwtPlot.xBottom, data.SpectrumE.min(), data.SpectrumE.max())

        plotData  = Qwt.QwtPlotCurve('Incident spectrum')
        plotData.setPen(Qt.QPen(Qt.Qt.red))
        
        plotData.setData(data.SpectrumE, data.Spectrum)
        plotData.attach(self.plot)
        
        layOut = QtGui.QVBoxLayout()
        layOut.addWidget(self.plot)

        self.setLayout(layOut)



class muPlot(QtGui.QWidget):
    def __init__(self, data, parent=None):
        QtGui.QWidget.__init__(self, parent)

        self.data = data

        self.vScaleSlider = Qwt.QwtSlider(self,Qt.Qt.Vertical) 
        self.vScaleSlider.setRange(data.muExt.min(), data.muExt.max(), 1.0)
        self.vScaleSlider.setValue(data.muExt.max())

        self.plotGrpB = QtGui.QGroupBox("Plot")
        self.plotGrpB.setMinimumHeight(250)
        self.plotGrpB.setMaximumWidth(200)

        self.showKAlpha  = QtGui.QCheckBox("K Alpha")
        self.showKBeta   = QtGui.QCheckBox("K Beta")
        self.showAbsEdge = QtGui.QCheckBox("Abs. edges")
        self.showIon     = QtGui.QCheckBox("Total Ion. coef.")
        self.showLine    = QtGui.QCheckBox("Line Ion. coef.")
        self.showRay     = QtGui.QCheckBox("Rayl. coef.")
        self.showTotal   = QtGui.QCheckBox("Total extinction")

        self.showKAlpha.setToolTip("Show K-alpha line markers")
        self.showKBeta.setToolTip("Show K-beta line markers")
        self.showAbsEdge.setToolTip("Show absorbtion edge markers")
        self.showIon.setToolTip("Show total photoionization coefficient")
        self.showLine.setToolTip("Show fluorescence line ionization coefficients")
        self.showRay.setToolTip("Show Rayleigh scattering coefficient")
        self.showTotal.setToolTip("Show total extinction coefficient")


        self.showKAlpha.setChecked(True)
        self.showAbsEdge.setChecked(True)
        self.showLine.setChecked(True)
        self.showIon.setChecked(True)
        self.showTotal.setChecked(True)

        plotTypeLayout = QtGui.QVBoxLayout()

        plotTypeLayout.addWidget(self.showKAlpha)
        plotTypeLayout.addWidget(self.showKBeta)
        plotTypeLayout.addWidget(self.showAbsEdge)
        plotTypeLayout.addWidget(self.showIon)
        plotTypeLayout.addWidget(self.showLine)
        plotTypeLayout.addWidget(self.showRay)
        plotTypeLayout.addWidget(self.showTotal)

        self.plotGrpB.setLayout(plotTypeLayout)

        self.muPlot = Qwt.QwtPlot()
        self.muPlot.setAxisScale(Qwt.QwtPlot.xBottom, data.energy.min(), data.energy.max())
        self.muPlot.setAxisScale(Qwt.QwtPlot.yLeft, 0, data.muExt.max())
        self.muPlot.setCanvasBackground(Qt.Qt.white)

        layOut = QtGui.QHBoxLayout()
        layOut.addWidget(self.muPlot)
        layOut.addWidget(self.vScaleSlider)
        layOut.addWidget(self.plotGrpB)
        self.setLayout(layOut)

        self.plot()

        self.connect(self.showKAlpha, Qt.SIGNAL('stateChanged(int)'), self.replot)
        self.connect(self.showKBeta, Qt.SIGNAL('stateChanged(int)'), self.replot)
        self.connect(self.showAbsEdge, Qt.SIGNAL('stateChanged(int)'), self.replot)
        self.connect(self.showIon, Qt.SIGNAL('stateChanged(int)'), self.replot)
        self.connect(self.showLine, Qt.SIGNAL('stateChanged(int)'), self.replot)
        self.connect(self.showRay, Qt.SIGNAL('stateChanged(int)'), self.replot)
        self.connect(self.showTotal, Qt.SIGNAL('stateChanged(int)'), self.replot)

        self.connect(self.vScaleSlider, Qt.SIGNAL('valueChanged(double)'), self.printVal)


    def plot(self):

        self.muExtCurve  = Qwt.QwtPlotCurve('Extinction coefficient')
        self.muExtCurve.setPen(Qt.QPen(Qt.Qt.black))
        self.muExtCurve.pen().setWidth(2)
        self.muExtCurve.setData(self.data.energy, self.data.muExt)
        self.muExtCurve.attach(self.muPlot)

        self.muPhotoIonCurve  = Qwt.QwtPlotCurve('Total Photoionization coefficient')
        self.muPhotoIonCurve.setPen(Qt.QPen(Qt.Qt.DashLine))
        self.muPhotoIonCurve.pen().setColor(Qt.Qt.red)
        self.muPhotoIonCurve.pen().setWidth(2)
        self.muPhotoIonCurve.setData(self.data.energy, self.data.muPhotoIon)
        self.muPhotoIonCurve.attach(self.muPlot)

        self.muIonCurve  = Qwt.QwtPlotCurve('Photoionization coefficient')
        self.muIonCurve.setPen(Qt.QPen(Qt.Qt.DashLine))
        self.muIonCurve.pen().setColor(Qt.Qt.black)
        self.muIonCurve.pen().setWidth(2)
        self.muIonCurve.setData(self.data.energy, self.data.muAbs.sum(1))
        self.muIonCurve.attach(self.muPlot)

        self.muRayCurve  = Qwt.QwtPlotCurve('Rayleigh scattering coefficient')
        self.muRayCurve.setPen(Qt.QPen(Qt.Qt.DashLine))
        self.muRayCurve.pen().setColor(Qt.Qt.blue)
        self.muRayCurve.pen().setWidth(4)
        self.muRayCurve.setData(self.data.energy, self.data.muRay)
        self.muRayCurve.attach(self.muPlot)

        self.muLineACurves = []
        self.muLineBCurves = []
        self.KaEdges      = []
        self.KaEdgeLabs   = []
        self.KbEdges      = []
        self.KbEdgeLabs   = []
        self.absEdges     = []

        for iElem in range(len(self.data.Elements)):

            self.muLineACurves.append(Qwt.QwtPlotCurve('Photoionization coefficient'))
            self.muLineACurves[iElem].setPen(Qt.QPen(Qt.Qt.DashLine))
            self.muLineACurves[iElem].pen().setColor(Qt.Qt.black)
            self.muLineACurves[iElem].pen().setWidth(1)
            self.muLineACurves[iElem].setData(self.data.energy, self.data.muAbs[:,2*iElem])
            self.muLineACurves[iElem].attach(self.muPlot)

            self.muLineBCurves.append(Qwt.QwtPlotCurve('Photoionization coefficient'))
            self.muLineBCurves[iElem].setPen(Qt.QPen(Qt.Qt.DashLine))
            self.muLineBCurves[iElem].pen().setColor(Qt.Qt.black)
            self.muLineBCurves[iElem].pen().setWidth(1)
            self.muLineBCurves[iElem].setData(self.data.energy, self.data.muAbs[:,2*iElem+1])
            self.muLineBCurves[iElem].attach(self.muPlot)

            self.KaEdges.append(Qwt.QwtPlotMarker())
            self.KaEdges[iElem].setLineStyle(Qwt.QwtPlotMarker.VLine)
            self.KaEdges[iElem].setLinePen(Qt.QPen(Qt.Qt.DashLine))
            self.KaEdges[iElem].linePen().setColor(Qt.Qt.red)
            self.KaEdges[iElem].setLabelAlignment(Qt.Qt.AlignRight | Qt.Qt.AlignTop)

            self.KaEdgeLabs.append(Qwt.QwtPlotMarker())
            self.KaEdgeLabs[iElem].setLineStyle(Qwt.QwtPlotMarker.NoLine)
            self.KaEdgeLabs[iElem].setLabelAlignment(Qt.Qt.AlignRight)
            self.KaEdgeLabs[iElem].setLabel(Qwt.QwtText(self.data.Elements[iElem] + ' Ka'))

            self.KaEdges[iElem].setXValue(self.data.lEnergy[iElem * 2])
            self.KaEdgeLabs[iElem].setValue(self.data.lEnergy[iElem * 2],  self.data.muExt.max() - 400*(iElem+1))

            self.KaEdges[iElem].attach(self.muPlot)
            self.KaEdgeLabs[iElem].attach(self.muPlot)


            self.KbEdges.append(Qwt.QwtPlotMarker())
            self.KbEdges[iElem].setLabelAlignment(Qt.Qt.AlignRight | Qt.Qt.AlignBottom)
            self.KbEdges[iElem].setLineStyle(Qwt.QwtPlotMarker.VLine)
            self.KbEdges[iElem].setLinePen(Qt.QPen(Qt.Qt.DotLine))
            self.KbEdges[iElem].linePen().setColor(Qt.Qt.blue)
            self.KbEdges[iElem].setXValue(self.data.lEnergy[iElem*2 + 1])
            self.KbEdges[iElem].attach(self.muPlot)

            self.KbEdgeLabs.append(Qwt.QwtPlotMarker())
            self.KbEdgeLabs[iElem].setLabelAlignment(Qt.Qt.AlignRight)
            self.KbEdgeLabs[iElem].setLineStyle(Qwt.QwtPlotMarker.NoLine)
            self.KbEdgeLabs[iElem].setLinePen(Qt.QPen(Qt.Qt.DotLine))
            self.KbEdgeLabs[iElem].linePen().setColor(Qt.Qt.blue)
            self.KbEdgeLabs[iElem].setLabel(Qwt.QwtText(self.data.Elements[iElem] + ' Kb'))
            self.KbEdgeLabs[iElem].setValue(self.data.lEnergy[iElem * 2 + 1],  400*(iElem+1))
            self.KbEdgeLabs[iElem].attach(self.muPlot)


            self.absEdges.append(Qwt.QwtPlotMarker())
            self.absEdges[iElem].setLabelAlignment(Qt.Qt.AlignRight | Qt.Qt.AlignCenter)
            self.absEdges[iElem].setLineStyle(Qwt.QwtPlotMarker.VLine)
            self.absEdges[iElem].setLinePen(Qt.QPen(Qt.Qt.SolidLine))
            self.absEdges[iElem].linePen().setColor(Qt.Qt.black)
            self.absEdges[iElem].setLabel(Qwt.QwtText(self.data.Elements[iElem]))
            self.absEdges[iElem].setXValue(self.data.eEnergy[iElem])
            self.absEdges[iElem].attach(self.muPlot)


        if(self.data.spcType == 'Line'):
            self.emissionLineEnergy = self.data.SpectrumE[0]

            self.emLine = Qwt.QwtPlotMarker()
            self.emLine.setLabelAlignment(Qt.Qt.AlignRight | Qt.Qt.AlignBottom)
            self.emLine.setLineStyle(Qwt.QwtPlotMarker.VLine)
            self.emLine.setLinePen(Qt.QPen(Qt.Qt.SolidLine))
            self.emLine.linePen().setColor(Qt.Qt.blue)
            self.emLine.setXValue(self.emissionLineEnergy)
            self.emLine.attach(self.muPlot)

        self.replot()

    def replot(self):

        self.muExtCurve.setVisible(self.showTotal.isChecked())
        self.muIonCurve.setVisible(self.showIon.isChecked())
        self.muRayCurve.setVisible(self.showRay.isChecked())
   
        for iElem in range(len(self.data.Elements)):
            self.muLineACurves[iElem].setVisible(self.showLine.isChecked())
            self.muLineBCurves[iElem].setVisible(self.showLine.isChecked())
            self.KaEdges[iElem].setVisible(self.showKAlpha.isChecked())
            self.KaEdgeLabs[iElem].setVisible(self.showKAlpha.isChecked())
            self.KbEdges[iElem].setVisible(self.showKBeta.isChecked())
            self.KbEdgeLabs[iElem].setVisible(self.showKBeta.isChecked())
            self.absEdges[iElem].setVisible(self.showAbsEdge.isChecked())

        self.muPlot.replot()

    def printVal(self, a):
        self.muPlot.setAxisScale(Qwt.QwtPlot.yLeft, 0, a)
        self.muPlot.replot()


class linePlot(QtGui.QWidget):
    def __init__(self, hs, data, parent=None):
        QtGui.QWidget.__init__(self, parent)

        self.data       = data
        self.ratioPlots = []
        self.theta      = np.linspace(0, 90, hs.resTheta)

        self.nElems     = len(hs.Elements)
        
        for i in range(self.nElems):

            self.ratioPlots.append( Qwt.QwtPlot())
            self.ratioPlots[i].enableAxis(Qwt.QwtPlot.xBottom, False)
            self.ratioPlots[i].setCanvasBackground(Qt.Qt.white)
            self.ratioPlots[i].setAxisScale(Qwt.QwtPlot.xBottom, 0, 90)
            #self.ratioPlots[i].setAxisScale(Qwt.QwtPlot.yLeft, self.dLims[i,0], self.dLims[i,1])
            self.ratioPlots[i].setAxisMaxMajor(Qwt.QwtPlot.yLeft, 4)
            self.ratioPlots[i].insertLegend(Qwt.QwtLegend(), Qwt.QwtPlot.RightLegend)

        self.plotData = []

        for iElem in range(self.nElems):
            self.plotData.append(Qwt.QwtPlotCurve(hs.Elements[iElem]))
            self.plotData[iElem].setPen(Qt.QPen(Qt.Qt.red))

            #I = self.data[0][iElem]
            I = self.data[0][iElem]

            self.plotData[iElem].setData(self.theta, I)
            self.plotData[iElem].attach(self.ratioPlots[iElem])

        layOut = QtGui.QVBoxLayout()

        for iElem in range(self.nElems):
            layOut.addWidget(self.ratioPlots[iElem])

        self.setMinimumWidth(400)
        self.setLayout(layOut)


    def replot(self, pType, iRefElem, iTht):
        
        for iElem in range(self.nElems):

            if(pType == 0):
                y = self.data[iTht][iElem]
                #self.ratioPlots[iElem].setAxisScale(Qwt.QwtPlot.yLeft, self.dLims[iElem,0], self.dLims[iElem,1])

            elif(pType == 1):
                self.ratioPlots[iElem].setAxisAutoScale(Qwt.QwtPlot.yLeft)
                y = self.data[iTht][iElem] / self.data[iTht][iRefElem]

            elif(pType == 2):
                self.ratioPlots[iElem].setAxisAutoScale(Qwt.QwtPlot.yLeft)
                y = self.data[iTht][iElem] / self.data[iTht][iRefElem]
                y /= y[0]

            self.plotData[iElem].setData(self.theta, y)
            self.plotData[iElem].attach(self.ratioPlots[iElem])

            self.ratioPlots[iElem].replot()  



class xrInspect(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)

        self.setWindowTitle('xrInspect')

        load_data(sys.argv[1], self)

        ## ELEMENT LIST
        ##
        self.elementGrpB = QtGui.QGroupBox("Elements")
        self.elementGrpB.setMinimumHeight(250)
        self.elementGrpB.setMaximumWidth(200)

        self.elementList = QtGui.QListWidget(self)
        self.elementList.addItems(self.hs.Elements)
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
        for i in range(self.hs.nThetaI):
            self.thetaList.addItem("%6.2f" % math.degrees(self.hs.thetaI[i]))
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
        self.linePlots = linePlot(self.hs, self.data, self)
        self.mu        = muPlot(self.hs, self)
        self.spectrum  = spectrumPlot(self.hs, self)
        self.hsPlots   = hsPlot(self.hs, self)

        self.hsPlots.setVisible(False)

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

        sbox.addWidget(self.linePlots)
        sbox.addWidget(self.hsPlots)

        
        infoLayout = QtGui.QVBoxLayout()
        infoLayout.addWidget(self.mu)
        if(self.hs.spcType=='Spectrum'):
            infoLayout.addWidget(self.spectrum)
        #infoLayout.addWidget(self.dPlot)

        sbox.addLayout(infoLayout)

        self.setLayout(sbox)

        self.resize(800, 1000)

        ## SET CONNECTIONS
        ##

        self.connect(self.plotStyleGroup, QtCore.SIGNAL('buttonClicked(int)'), self.plot)

        self.connect(self.elementList, QtCore.SIGNAL('itemSelectionChanged()'), self.plot)
        self.connect(self.thetaList, QtCore.SIGNAL('itemSelectionChanged()'), self.plot)

        self.connect(self.showPlots, QtCore.SIGNAL('stateChanged(int)'), self.togglePlots)
        self.connect(self.showInfo, QtCore.SIGNAL('clicked()'), self.toggleInfo)
        self.connect(self.qb, QtCore.SIGNAL('clicked()'), QtGui.qApp, QtCore.SLOT('quit()'))

        self.connect(self.hb, QtCore.SIGNAL('clicked()'), self.hide)

    def plot(self):

        pType    = self.plotStyleGroup.checkedId()
        iRefElem = self.elementList.currentRow()
        iTht     = self.thetaList.currentRow()

        if(self.linePlots.isVisible()):
            self.linePlots.replot(pType, iRefElem, iTht)

        if(self.hsPlots.isVisible()):
            self.hsPlots.replot(pType, iRefElem, iTht)

        if pType > 0:
            self.elementList.setEnabled(True)
        else:
            self.elementList.setEnabled(False)

    def togglePlots(self, v):
        self.linePlots.setVisible(v)

    def toggleInfo(self):
        self.spectrum.setVisible(self.showInfo.isChecked())
        self.mu.setVisible(self.showInfo.isChecked())


    def hide(self):
        if self.linePlots.isVisible():
            self.linePlots.setVisible(False)
            self.hsPlots.setVisible(True)
            self.plot()
        else:
            self.linePlots.setVisible(True)
            self.hsPlots.setVisible(False)
            self.plot()


def load_data(fileName, w):
    
    w.hs = xrHemisphere() 
    w.hs.load(fileName)

    if(w.hs.Method == "MonteCarlo"):
        w.hs.divideBySolidAngle()
    
    w.data          = []
    w.dLims         = np.zeros([len(w.hs.Elements), 2])

    for iTht in range(w.hs.nThetaI):
        w.data.append([])
        for iElem in range(len(w.hs.Elements)):
            w.data[iTht].append((w.hs.toArray(set=iTht, lvl=iElem*2) + w.hs.toArray(set=iTht, lvl=iElem*2+1)).mean(1))

    if 'Analytic' in w.hs.hsData:
        w.data_analytic = []
        
        for iTht in range(w.hs.nThetaI):
            w.data_analytic.append([])
            for iElem in range(len(w.hs.Elements)):
                w.data_analytic[iTht].append((w.hs.toArray(set=iTht, lvl=iElem*2, hsDataSet='Analytic') + w.hs.toArray(set=iTht, lvl=iElem*2+1, hsDataSet='Analytic')).mean(1))

    for iElem in range(len(w.hs.Elements)):
        minV = 1e8
        maxV = -1e8

        for iTht in range(w.hs.nThetaI):
            if w.data[iTht][iElem].min() < minV:
                minV = w.data[iTht][iElem].min()
            if w.data[iTht][iElem].max() > maxV:
                maxV = w.data[iTht][iElem].max()

        w.dLims[iElem, :] = [minV, maxV]


app = QtGui.QApplication(sys.argv)
win = xrInspect()

win.show()
app.exec_()
