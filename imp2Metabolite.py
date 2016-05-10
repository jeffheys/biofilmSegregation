# -*- coding: utf-8 -*-
"""
Created on Tue Dec 30 09:15:59 2014

@author: jheys

Notes:
    Domain size is in millimeters. Should be ~0.2 in order to be ~200 times bug radius
    Mesh size should be about 20x20
    Diffusivity is in mm^2/h to be approximately order zero for most substrates
"""

import numpy
import matplotlib.pyplot as plt
import scipy.linalg as slin

class substrate:
    def __init__(self, name, sizeDomain, meshSize, diff, bulk, dt):
        self.name = name
        self.sizeDomain = sizeDomain
        self.meshSize = meshSize
        self.xLoc = numpy.linspace(0.0,self.sizeDomain[0],num=self.meshSize[0])
        self.dx = self.xLoc[1] - self.xLoc[0]
        self.yLoc = numpy.linspace(0.0,self.sizeDomain[1],num=self.meshSize[1])
        self.dy = self.yLoc[1] - self.yLoc[0]
        self.pts = meshSize[0]*meshSize[1]
        self.conc = numpy.zeros((self.pts,),dtype=numpy.float)
        self.diff = -0.5 * diff * dt
        self.bulk = bulk
        for i in range(self.meshSize[0]):  # intitial conditions on top edge
            node = self.nodeNum(i,self.meshSize[1]-1)
            self.conc[node] = self.bulk
        # Build diffusion operator
        self.opA = numpy.zeros((self.pts,self.pts),dtype=numpy.float)
        for j in range(self.meshSize[1]): # loop through y nodes
            for i in range(self.meshSize[0]): # loop through x nodes
                node = self.nodeNum(i,j)
                nodeS = self.nodeNum(i,j-1)
                nodeE = self.nodeNum(i+1,j)
                nodeN = self.nodeNum(i,j+1)
                nodeW = self.nodeNum(i-1,j)
                bcX = False
                bcY = False
                self.opA[node,node] += 1.0
                if i == 0: # west wall
                    bcX = True
                    self.opA[node,node] += -2.0*self.diff/(self.dx**2)
                    self.opA[node,nodeE] += 2.0*self.diff/(self.dx**2)
                if i == (self.meshSize[0]-1): # east wall
                    bcX = True
                    self.opA[node,node] += -2.0*self.diff/(self.dx**2)
                    self.opA[node,nodeW] += 2.0*self.diff/(self.dx**2)
                if j == 0: # south wall
                    bcY = True
                    self.opA[node,node] += -2.0*self.diff/(self.dy**2)
                    self.opA[node,nodeN] += 2.0*self.diff/(self.dy**2)
                if j == (self.meshSize[1]-1): # top surface
                    bcX = True
                    bcY = True
                    self.opA[node,node] = 1.0
                    if i == 0:
                        self.opA[node,nodeE] = 0.0
                    if i == (self.meshSize[0]-1):
                        self.opA[node,nodeW] = 0.0
                if bcX == False:
                    self.opA[node,node] += -2.0*self.diff/(self.dx**2)
                    self.opA[node,nodeW] += 1.0*self.diff/(self.dx**2)
                    self.opA[node,nodeE] += 1.0*self.diff/(self.dx**2)
                if bcY == False:
                    self.opA[node,node] += -2.0*self.diff/(self.dy**2)
                    self.opA[node,nodeS] += 1.0*self.diff/(self.dy**2)
                    self.opA[node,nodeN] += 1.0*self.diff/(self.dy**2)
        self.lu = slin.lu_factor(self.opA)

    def nodeNum(self, i, j):
        return j*self.meshSize[0] + i

    def diffuse(self):  # uses Crank-Nicolson time stepping
        b = numpy.zeros(self.pts,dtype=numpy.float)
        for j in range(self.meshSize[1]): # loop through y nodes
            for i in range(self.meshSize[0]): # loop through x nodes
                node = self.nodeNum(i,j)
                nodeS = self.nodeNum(i,j-1)
                nodeE = self.nodeNum(i+1,j)
                nodeN = self.nodeNum(i,j+1)
                nodeW = self.nodeNum(i-1,j)
                bcX = False
                bcY = False
                b[node] += self.conc[node]
                if i == 0: # west wall
                    bcX = True
                    b[node] += 2.0*self.diff*self.conc[node]/(self.dx**2)
                    b[node] += -2.0*self.diff*self.conc[nodeE]/(self.dx**2)
                if i == (self.meshSize[0]-1): # east wall
                    bcX = True
                    b[node] += 2.0*self.diff*self.conc[node]/(self.dx**2)
                    b[node] += -2.0*self.diff*self.conc[nodeW]/(self.dx**2)
                if j == 0: # south wall
                    bcY = True
                    b[node] += 2.0*self.diff*self.conc[node]/(self.dy**2)
                    b[node] += -2.0*self.diff*self.conc[nodeN]/(self.dy**2)
                if j == (self.meshSize[1]-1): # top surface
                    bcX = True
                    bcY = True
                    b[node] = self.bulk
                if bcX == False:
                    b[node] += 2.0*self.diff*self.conc[node]/(self.dx**2)
                    b[node] += -1.0*self.diff*self.conc[nodeW]/(self.dx**2)
                    b[node] += -1.0*self.diff*self.conc[nodeE]/(self.dx**2)
                if bcY == False:
                    b[node] += 2.0*self.diff*self.conc[node]/(self.dy**2)
                    b[node] += -1.0*self.diff*self.conc[nodeS]/(self.dy**2)
                    b[node] += -1.0*self.diff*self.conc[nodeN]/(self.dy**2)
        self.conc = slin.lu_solve(self.lu,b)
        return True

    def showPlot(self,figNum=0,bugs=0):
        ax = plt.figure(figNum)
        plt.xlim((0.0,self.sizeDomain[0]))
        plt.ylim((0.0,self.sizeDomain[1]))
        X,Y = numpy.meshgrid(self.xLoc,self.yLoc)
        maxC = "Maximum = %.4f" % numpy.amax(self.conc)
        avgC = "Average = %.4f" % numpy.average(self.conc)
        numBugs = 'Number of bugs = %s' % bugs
        for j in range(self.meshSize[1]): # loop through nodes in the y-direction
            for i in range(self.meshSize[0]):
                node = self.nodeNum(i,j)
                if self.conc[node] < 0.0:
                    print "Warning: negative", self.name, "concentration at node ", node
                    print 'i = ', i, ', j = ', j
                    self.conc[node] = 0.0
        cont1=plt.contour(X,Y,self.conc.reshape((self.meshSize[0],self.meshSize[1])))
        cb=plt.colorbar(cont1)
        cb.set_label(self.name.title())
        ax.text(0.5, 0.82, maxC)
        ax.text(0.5, 0.85, avgC)
        ax.text(0.5, 0.79, numBugs)
        for line in cb.lines:
            line.set_linewidth(10)


if __name__ == "__main__":
    import time
    start = time.time()
    tstep = 0.1/60.0
    ace = substrate('acetate', sizeDomain = (0.2,0.2), meshSize = (24,24),
diff=10., bulk=1., dt=tstep)
    for t in range(50):
        ace.diffuse()
    ace.showPlot(figNum=1)
    print 'End concentrations:\n',ace.conc
    end = time.time()
    print 'Run time:', (end-start)/60, ' minutes'
