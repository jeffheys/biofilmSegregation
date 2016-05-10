# -*- coding: utf-8 -*-
"""
Created on Tue Dec 30 09:17:16 2014

@author: jheys
Edited by Michael Venters

Notes:
    Mass in units of picograms (cell dry weight)
    .28 pg corresponds to about 1 cubic micron cell volume (rule of thumb)
"""

class Ecoli:
    def __init__(self, location, pheno='pp', radius=.001, mass=.28, divideFactor=2.0, orig = 0):
        self.location = location
        self.pheno = pheno
        self.orig = orig
        self.radius = radius
        self.mass = mass
        self.initialMass = mass
        self.divideMass = mass * divideFactor
        if self.pheno == 'pp':  #primary producer
            self.mu = 0.6 # per hour
            self.Kg = 0.05 # g/L
            self.Kip = 0.25 # g/L
            self.Y11 = 0.152 # mg cdw/mg glucose consumed
            self.Y21 = 0.2 # mg cdw/mg acetate produced
        elif self.pheno == 'fixed':   #producer with fixed growth rate
            self.mu = 0.6 # per hour
            self.Kg = 0.05 # g/L
            self.Kip = 0.5 # g/L
            self.Y11 = 0.152 # mg cdw/mg glucose consumed
            self.Y21 = 0.114 # mg cdw/mg acetate produced
        elif self.pheno == 's1':  #scavenger
            self.mu = 0.3 # per hour
#            self.Ka = 0.23 # g/L
#            self.Kis = 4.6 # g/L
            self.Ka = 0.05 # g/L
            self.Kis = 0.7 # g/L
            self.KO2 = 2.75e-4 # g/L
            self.Y22 = 0.40 # mg cdw/mg acetate consumed
            self.YO2 = 1.0 # mg cdw/mg O2 consumed

    def grow(self, G, A, O2, dt):
        if self.pheno == 'pp':
            growthRate = self.mu * G / (G + self.Kg) # glucose effects, h^-1
            growthRate = growthRate * self.Kip / (self.Kip + A)  # inhibition
#            print 'Producer specific growth rate:', growthRate, ' h^-1'
#            print 'Glucose conc:', G, ' g/L'
            self.mass += growthRate * self.mass * dt
        elif self.pheno == 'fixed':   # fixed growth rate, no size change
            growthRate = 0.6
        elif self.pheno == 's1':
            growthRate = self.mu*A*O2 / ((A + self.Ka + A * A / self.Kis)*(self.KO2 + O2))  # h^-1
#            print 'Scavenger specific growth rate:', growthRate, ' h^-1'
#            print 'Acetate conc:', A, ' g/L'
#            print 'O2 conc:', O2, ' g/L'
            self.mass += growthRate * self.mass * dt
        growthRate = growthRate*self.mass    # dx/dt [=] picogram cdw/h
#        print 'dx/dt =', growthRate, ' pg cdw/h'
        return growthRate

    def divide(self):
        return self.mass > self.divideMass


if __name__ == '__main__':
    bug = Ecoli([.5,.5],pheno='s1')
    print(bug.grow(2.4,.47,0.1,1/60.))
