# -*- coding: utf-8 -*-
"""
Created on Tue Dec 30 09:15:59 2014

@author: jheys
Edited by Michael Venters

Notes:
    Domain size in millimeters
    Diffusivities in mm^2/h
    K values and concentrations in g/L
"""

from bugs import Ecoli
from impMetabolite import substrate
import numpy
import matplotlib.pyplot as plt
import random
from time import time, asctime
import sys

class Biofilm:
    def __init__(self, sizeDomain = (0.2,0.2), meshSize = (32,32), numS1 = 16, clusterSize=0.1, diffDt = 0.1, Adiff=0.994):
        self.sizeDomain = sizeDomain
        self.areaScale = numpy.power(sizeDomain[0] / 0.2, 1.5)
        self.meshSize = meshSize
        self.xLoc = numpy.linspace(0.0,self.sizeDomain[0],num=self.meshSize[0])
        self.dx = self.xLoc[1] - self.xLoc[0]
        self.yLoc = numpy.linspace(0.0,self.sizeDomain[1],num=self.meshSize[1])
        self.dy = self.yLoc[1] - self.yLoc[0]
        self.vcube = self.dx * self.dy * self.dx # millimeter^3

        #make initial population of bugs
        self.bugList = []
        radius = 0.001   # mm (1 micron)
        mass = 1.5*0.28  # picograms  (median mass to division)
        pheno = 'fixed'    #make producer in the center of the domain
        x = sizeDomain[0]*0.5
        y = sizeDomain[1]*0.5
        bug = Ecoli([x,y],pheno,radius,mass)
        self.bugList.append(bug)
        for i in range(numS1):    #randomly cluster numS1 scavenger bugs around the producer
            radius = 0.001   # mm (1 micron)
            mass = 0.28  # picograms
            pheno = 's1'
            r = clusterSize * sizeDomain[0] * random.random() # clusterSize is a percentage of domain size
            theta = 2.0 * numpy.pi * random.random()
            x = sizeDomain[0]*0.5 + r * numpy.cos(theta)
            y = sizeDomain[1]*0.5 + r * numpy.sin(theta)
            bug = Ecoli([x,y],pheno,radius,mass,orig=i)      # NOTE: first scavenger will have orig=0
            self.bugList.append(bug)
        self.foods = []
        newfood = substrate('glucose', self.sizeDomain, self.meshSize,diff=0.84,bulk=1.0, dt=diffDt)
        self.foods.append(newfood)
        newfood = substrate('acetate', self.sizeDomain, self.meshSize,diff=Adiff,bulk=0.0, dt=diffDt)
        self.foods.append(newfood)
        newfood = substrate('oxygen', self.sizeDomain, self.meshSize,diff=9.36,bulk=8.32e-3, dt=diffDt)
        self.foods.append(newfood)

    def spread(self):
        numBugs = len(self.bugList)
        numMoves = 1
        maxItr = 10
        itr=0
        while numMoves > 0 and itr < maxItr:
            itr += 1
            numMoves = 0
            for i in range(numBugs):
                for j in range(numBugs):
                    if i != j:
                        dx = (self.bugList[i].location[0] - self.bugList[j].location[0])
                        dy = (self.bugList[i].location[1] - self.bugList[j].location[1])
                        dist = numpy.sqrt(dx**2 + dy**2)
                        if dist < 1.0e-14:
                            dist = 1.0e-14
                            print("ERROR: two bugs at the same location!")
                        allowDist = self.bugList[i].radius + self.bugList[j].radius
                        if dist < allowDist:   # check for collision
                            numMoves += 1
                            dirMove = (dx/dist, dy/dist)
                            moveSize = (allowDist - dist) + random.random() * allowDist
                            self.bugList[i].location[0] += dirMove[0]*moveSize
                            self.bugList[i].location[1] += dirMove[1]*moveSize
                            #check right boundary
                            if (self.bugList[i].location[0] + self.bugList[i].radius) > self.sizeDomain[0]:
                                self.bugList[i].location[0] = self.sizeDomain[0] - self.bugList[i].radius
                            #check left boundary
                            if (self.bugList[i].location[0] - self.bugList[i].radius) < 0.0:
                                self.bugList[i].location[0] = self.bugList[i].radius
                            #check top boundary
                            if (self.bugList[i].location[1] + self.bugList[i].radius) > self.sizeDomain[1]:
                                print('ERROR: Bugs above the top of the domain!')
                            #check bottom boundary
                            if (self.bugList[i].location[1] - self.bugList[i].radius) < 0.0:
                                self.bugList[i].location[1] = self.bugList[i].radius
#            print 'Iteration: ', itr, ' with ', numMoves, ' moves'
#            self.showPlot(figNum=itr)
#        print('spread complete in ', itr, ' iterations')
        
    def eat(self, dt):
        for b in self.bugList:
            SW = (int(b.location[0]/self.dx),int(b.location[1]/self.dy))
            nodeSW = SW[1]*self.meshSize[0] + SW[0]
            nodeSE = SW[1]*self.meshSize[0] + SW[0] + 1
            nodeNW = (SW[1]+1)*self.meshSize[0] + SW[0]
            nodeNE = (SW[1]+1)*self.meshSize[0] + SW[0] + 1
            dx1 =  b.location[0] - self.xLoc[SW[0]]
            dy1 = b.location[1] - self.yLoc[SW[1]]
            dx2 = self.dx - dx1
            dy2 = self.dy - dy1
            distSW = dx2*dy2/(self.dx*self.dy)
            distSE = dx1*dy2/(self.dx*self.dy)
            distNW = dx2*dy1/(self.dx*self.dy)
            distNE = dx1*dy1/(self.dx*self.dy)
#            print distSW, distSE, distNW, distNE
            for food in self.foods:
                cSW = food.conc[nodeSW]
                cSE = food.conc[nodeSE]
                cNW = food.conc[nodeNW]
                cNE = food.conc[nodeNE]
                #approx. substrate conc. at bug's location
                conc = cSW * distSW + cSE * distSE + cNW * distNW + cNE * distNE
                if conc < 0.0:
                    print('***Warning:', food.name,' concentration < 0!')
#                print cSW, cSE, cNW, cNE
                if food.name == 'glucose':
                    G = conc
                if food.name == 'acetate':
                    A = conc
                if food.name == 'oxygen':
                    O2 = conc
            growthRate = b.grow(G,A,O2,dt)  # returns units of pg cdw/h ~order 10^-1
            if b.pheno == 'fixed':
                for food in self.foods:
                    eaten = 0.0
                    if food.name == 'glucose':
                        dG = growthRate / b.Y11 /self.vcube * dt  # pg/mm^3 of glucose consumed
                        dG = dG * 10**-6   # convert to g/L
#                        print 'dG:',-dG
                        if food.conc[nodeSW] > dG * distSW:
                            food.conc[nodeSW] -= dG * distSW
                            eaten += dG * distSW
                        else:
                            eaten += food.conc[nodeSW]
                            food.conc[nodeSW] = 0.0
                        if food.conc[nodeSE] > dG * distSE:
                            food.conc[nodeSE] -= dG * distSE
                            eaten += dG * distSE
                        else:
                            eaten += food.conc[nodeSE]
                            food.conc[nodeSE] = 0.0
                        if food.conc[nodeNW] > dG * distNW:
                            food.conc[nodeNW] -= dG * distNW
                            eaten += dG * distNW
                        else:
                            eaten += food.conc[nodeNW]
                            food.conc[nodeNW] = 0.0
                        if food.conc[nodeNE] > dG * distNE:
                            food.conc[nodeNE] -= dG * distNE
                            eaten += dG * distNE
                        else:
                            eaten += food.conc[nodeNE]
                            food.conc[nodeNE] = 0.0
#                        if eaten < 0.999*dG:
#                            print "     Warning: Glucose is limiting (", eaten, "eaten)"
                    if food.name == 'acetate':
                        dA = growthRate / b.Y21 /self.vcube * dt * 10**-6  #convert pg/mm^3 to g/L
                        dA = dA * self.areaScale
#                        print 'fixed dA:', dA
                        food.conc[nodeSW] += dA * distSW
                        food.conc[nodeSE] += dA * distSE
                        food.conc[nodeNW] += dA * distNW
                        food.conc[nodeNE] += dA * distNE
            if b.pheno == 's1':
                for food in self.foods:
                    if food.name == 'oxygen':
                        eaten = 0.0
                        dO = growthRate / b.YO2 /self.vcube * dt * 10**-6  #convert pg/mm^3 to g/L
                        if food.conc[nodeSW] > dO * distSW:
                            food.conc[nodeSW] -= dO * distSW
                            eaten += dO * distSW
                        else:
                            eaten += food.conc[nodeSW]
                            food.conc[nodeSW] = 0.0
                        if food.conc[nodeSE] > dO * distSE:
                            food.conc[nodeSE] -= dO * distSE
                            eaten += dO * distSE
                        else:
                            eaten += food.conc[nodeSE]
                            food.conc[nodeSE] = 0.0
                        if food.conc[nodeNW] > dO * distNW:
                            food.conc[nodeNW] -= dO * distNW
                            eaten += dO * distNW                            
                        else:
                            eaten += food.conc[nodeNW]
                            food.conc[nodeNW] = 0.0
                        if food.conc[nodeNE] > dO * distNE:
                            food.conc[nodeNE] -= dO * distNE
                            eaten += dO * distNE
                        else:
                            eaten += food.conc[nodeNE]
                            food.conc[nodeNE] = 0.0
#                        if eaten < 0.999*dO:
#                            print "    Warning: Oxygen is limiting (", eaten, "eaten)"
                    if food.name == 'acetate':
                        eaten = 0.0
                        dA = growthRate / b.Y22 /self.vcube * dt * 10**-6  #convert pg/mm^3 to g/L
#                        print 's1 dA:', -dA
                        if food.conc[nodeSW] > dA * distSW:
                            food.conc[nodeSW] -= dA * distSW
                            eaten += dA * distSW
                        else:
                            eaten += food.conc[nodeSW]
                            food.conc[nodeSW] = 0.0
                        if food.conc[nodeSE] > dA * distSE:
                            food.conc[nodeSE] -= dA * distSE
                            eaten += dA * distSE
                        else:
                            eaten += food.conc[nodeSE]
                            food.conc[nodeSE] = 0.0
                        if food.conc[nodeNW] > dA * distNW:
                            food.conc[nodeNW] -= dA * distNW
                            eaten += dA * distNW                            
                        else:
                            eaten += food.conc[nodeNW]
                            food.conc[nodeNW] = 0.0
                        if food.conc[nodeNE] > dA * distNE:
                            food.conc[nodeNE] -= dA * distNE
                            eaten += dA * distNE
                        else:
                            eaten += food.conc[nodeNE]
                            food.conc[nodeNE] = 0.0
#                        if eaten < 0.999*dA:
#                            print "    Warning: Acetate is limiting (", eaten, "eaten)"

    def clone(self):
        for b in self.bugList:
            if b.divide() == True:
#                print('dividing')
                b.mass = b.initialMass
                x = b.location[0] + random.choice([-1,1])*(1.0 + random.random())*b.radius*2.0
                y = b.location[1] + random.choice([-1,1])*(1.0 + random.random())*b.radius*2.0
                bug = Ecoli([x,y],b.pheno,b.radius,b.mass)
                self.bugList.append(bug)

    def timeLoop(self, timeSteps, diffTimeSteps, growDt):
        for t in range(timeSteps):
            print("Time Step: ", t)
            for food in self.foods:
                try:
                    for i in range(diffTimeSteps):
                        food.diffuse()
                except:
                    print('***A diffusion error occured')
                    print(sys.exc_info())
                    return -1.
            try:
                self.eat(growDt)
            except IndexError as e:
                print('\n**Bugs above the top of the domain. Stopping')
                return t
            if t%10 == 0:
                random.shuffle(self.bugList)
            self.clone()
            self.spread()
#            if t%60 == 0:
#                self.showPlot(figNum=t,show=['acetate'],bugs=len(self.bugList))
            print('Total Bugs: ', len(self.bugList))
        return len(self.bugList)

    def showPlot(self,figNum=0,show=[],bugs=0):
        fig = plt.figure(figNum)
        plt.xlim((0.0,self.sizeDomain[0]))
        plt.ylim((0.0,self.sizeDomain[1]))
        scale = 20000.0
        bugColors = {'pp':'b', 's1':'r', 'fixed':'g'}
        for b in self.bugList:
            plt.scatter(b.location[0],b.location[1],s=b.radius*scale,c=bugColors[b.pheno])
        for name in show:
            for food in self.foods:
                if name == food.name:
                    food.showPlot(figNum,bugs)
        filename = 'singleCluster' + str(figNum) + '.png'
        fig.savefig(filename)

    def segregation(self):
        numMatch = 5
        totalMatch = 0
        totalCompare = 0
        for b in self.bugList:
            matches = []
            for i in range(numMatch):
                matches.append([1.0, 'xx'])
            for b2 in self.bugList:
                dx = b.location[0] - b2.location[0]
                dy = b.location[1] - b2.location[1]
                dist = numpy.sqrt(dx**2 + dy**2)
                for i in range(numMatch):
                    if dist < matches[i][0] and dist > 0.0:
                        matches[i][0] = dist
                        matches[i][1] = b2.pheno
                        break
            for i in range(numMatch):
                if b.pheno == matches[i][1]:
                    totalMatch += 1
                totalCompare += 1
        print(b.pheno, " at ",b.location, matches)
        return (totalMatch*1.0)/totalCompare

start = time()
numRuns = 6
numStartCells = 16    #number of initial scavengers (S1)
endtime = 150.0  # hours
timesteps = int(endtime*30) # < 1-minute grow steps
diffTimeSteps = 50 # number of diffusion steps per growth step
growDt = 1.0*endtime / timesteps
diffDt = 1.0*endtime / (timesteps*diffTimeSteps)

runSizes = [0.05, 0.15, 0.25]    #cluster sizes as a percentage of domain size
runDiffs = [0.2, 0.5, 2.0, 5.0, 20.0, 50.0, 200.0] # domain size
grow = numpy.zeros((len(runSizes)+1, len(runDiffs)+1))
growStd = numpy.zeros((len(runSizes)+1, len(runDiffs)+1))
maxDiv = numpy.zeros((len(runSizes)+1, len(runDiffs)+1))
minDiv = numpy.zeros((len(runSizes)+1, len(runDiffs)+1))
sizeNum = 0
for size in runSizes:
    diffNum = 0
    grow[sizeNum+1,0]=size
    growStd[sizeNum+1,0]=size
    maxDiv[sizeNum+1,0]=size
    minDiv[sizeNum+1,0]=size
    for diff in runDiffs:
        grow[0,diffNum+1]=diff
        growStd[0,diffNum+1]=diff
        maxDiv[0,diffNum+1]=diff
        minDiv[0,diffNum+1]=diff
        runGrow = numpy.zeros(numRuns)
        runMaxDiv = numpy.zeros(numRuns)
        runMinDiv = numpy.zeros(numRuns)
        for run in range(numRuns):
            parent = numpy.zeros(numStartCells)
            biofilmSim = Biofilm(sizeDomain=(diff,diff),numS1=numStartCells, clusterSize=size, diffDt=diffDt)
            biofilmSim.spread()
            runGrow[run] = biofilmSim.timeLoop(timesteps, diffTimeSteps, growDt)
            for b in biofilmSim.bugList:
                if b.pheno != 'fixed':
                    parent[b.orig] += 1    #count # of descendents from each original bug (scavengers only)
            runMaxDiv[run] = numpy.amax(parent)
            runMinDiv[run] = numpy.amin(parent)
        grow[sizeNum+1,diffNum+1]=numpy.average(runGrow)
        growStd[sizeNum+1,diffNum+1]=numpy.std(runGrow)
        maxDiv[sizeNum+1,diffNum+1]=numpy.average(runMaxDiv)
        minDiv[sizeNum+1,diffNum+1]=numpy.average(runMinDiv)
        biofilmSim.showPlot(figNum=(sizeNum*10+diffNum),show=['acetate'], bugs=len(biofilmSim.bugList))
        diffNum += 1
    sizeNum += 1

numpy.savetxt("grow.csv", grow, delimiter=",")
numpy.savetxt("growStd.csv", growStd, delimiter=",")
numpy.savetxt("maxDiv.csv", maxDiv, delimiter=",")
numpy.savetxt("minDiv.csv", minDiv, delimiter=",")

end = time()
print("Simulation time:", endtime, "hours")
print('Elapsed time:', (end-start)/60, ' minutes')
print('Exiting.')
print(asctime())
