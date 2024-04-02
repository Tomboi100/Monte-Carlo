# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 14:47:48 2017

@author: edf
"""

import numpy as np
import numpy.random as rnd
import matplotlib.pyplot as plt

def RandomVector(d):
    return 2*rnd.random(d)-1

def IsotropicVector(d):
    while True:
        u = RandomVector(d)
        if(np.sum(u**2)<1):
            break
    return u

def AllocateAtoms(Natoms,d):
    return np.empty([Natoms,d])

def InitPositions(r,Natoms,d,boxsize):
    for i in range(Natoms):
        r[i,:] = IsotropicVector(d)*boxsize

def InsideSphere(p,d,boxsize):
    return np.sum(p**2)<boxsize**2

def AllInsideSphere(r,Natoms,d,boxsize):
    inside = True
    for i in range(Natoms):
        inside = inside and InsideSphere(r[i,:],d,boxsize)
    return inside

def GenerateMove(r,Natoms,d,stepsize):
    imov = rnd.randint(Natoms)
    oldpos = r[imov,:]
    u = IsotropicVector(d)
    newpos = oldpos + stepsize*u
    return imov, oldpos, newpos

def AcceptMove(imov,oldpos,newpos,r,Natoms,d):
    r[imov,:] = newpos

def RejectMove(imov,oldpos,newpos,r,Natoms,d):
    r[imov,:] = oldpos

def Potential(r2):
    return 4*((1.0/r2**6)-(1.0/r2**3))

def ComputeEnergy(r,Natoms,d):
    E=0.0
    for i in range(Natoms-1):
        for j in range(i+1,Natoms):
            r2 = np.sum((r[i,:]-r[j,:])**2)
            E += Potential(r2)
    return E

def ComputeDiffEnergy(imov,oldpos,newpos,r,Natoms,d):
    DE = 0.0
    for i in range(Natoms):
        if (not(i == imov)):
            r2 = np.sum((r[i,:]-newpos)**2)
            DE += Potential(r2)
            r2 = np.sum((r[i,:]-oldpos)**2)
            DE -= Potential(r2)
    return DE

def MCstep(r,Natoms,d,Energy,nacc,Temperature,boxsize,stepsize):
    imov, oldpos, newpos = GenerateMove(r,Natoms,d,stepsize)
    if (not InsideSphere(newpos,d,boxsize)):
        RejectMove(imov,oldpos,newpos,r,Natoms,d)
    else:
        DE = ComputeDiffEnergy(imov,oldpos,newpos,r,Natoms,d)
        if DE <= 0.0:
            AcceptMove(imov,oldpos,newpos,r,Natoms,d)
            Energy += DE
            nacc += 1
        else:
            if rnd.random()<np.exp(-DE/Temperature):
                AcceptMove(imov,oldpos,newpos,r,Natoms,d)
                Energy += DE
                nacc += 1
            else:
                RejectMove(imov,oldpos,newpos,r,Natoms,d)
    return Energy, nacc

def PlotPositions(r,Natoms,d,boxsize):
    plt.scatter(r[:,0],r[:,1])
    plt.show()

d=2
boxsize = 2.5
Temperature = 1.0
stepsize = 0.5
Natoms = 10
Niter = 100000
Nequil = Niter//10
r = AllocateAtoms(Natoms,d)
InitPositions(r,Natoms,d,boxsize)

# Equilibration:
nacc = 0
Energy = ComputeEnergy(r,Natoms,d)
for i in range(Nequil):
    Energy, nacc = MCstep(r,Natoms,d,Energy,nacc,
                          Temperature,boxsize,stepsize)

# Production run
nacc = 0
sum_Energy = 0.0
sum_Energy2 = 0.0
Energy = ComputeEnergy(r,Natoms,d)
for i in range(Niter):
    Energy, nacc = MCstep(r,Natoms,d,Energy,nacc,
                          Temperature,boxsize,stepsize)
    sum_Energy += Energy
    sum_Energy2 += Energy**2

# Calculate averages
av_Energy = sum_Energy/float(Niter)
av_Energy2 = sum_Energy2/float(Niter)
CV = (av_Energy2-av_Energy**2)/Temperature**2

# print results
print("acceptance ratio",nacc/float(Niter))
print("average energy",av_Energy)
print("heat capacity",CV)

PlotPositions(r,Natoms,d,boxsize)
