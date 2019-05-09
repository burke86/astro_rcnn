#!/usr/bin/env python

"""
  @package phosim
  @file prescriptionDataToPhosim.py
  @brief python script to convert a zemax prescription data report file to phosim perturbation.txt body commands.
  This script will correctly convert non-axial optical systems (with rotations and translations) using global
  coordinate system data from Zemax (Analyze->Reports->Prescription Data). 
 
  @brief Created by:
  @author Colin Burke (Purdue)
 
  @warning This code is not fully validated
  and not ready for full release.  Please
  treat results with caution.

Usage: python prescriptionDataToPhosim.py config.txt zemaxPDR.txt perturbation.txt

       config.txt is a user provided file which contains
       Column 1: surface name (primary, secondary, etc.) 
       Column 2: Zemax surface number (from prescription data file) 
       Column 3: PhoSim surface number (surface link from optics file) 

       e.g.,

       primary    4 0
       secondary  5 1
       tertiary   9 2


Notes: 1. This script assumes data in mm and degrees, and the global reference surface lies at the origin and is not rotated. 
       2. It will overwrite existing perturbation.txt in the current folder.
       3. requires external module 'transformations.py' by Christoph Gohlke, <http://www.lfd.uci.edu/~gohlke/code/transformations.py.html>.

"""

import sys, os, commands
import math as m
import transformations as tf


class BodyCommand(object):
    def  __init__(self,dx,dy,dz,rx,ry,rz,surfID):
        self.dx=dx
        self.dy=dy
        self.dz=-dz
        self.surfID=surfID #surface id corresponding to names, zmxSurfNos, and phosimSurfNos lists.
        M=tf.euler_matrix(m.radians(rx),m.radians(ry),m.radians(rz),'rxyz')
        phi,theta,psi=tf.euler_from_matrix(M,'rzxz')
        theta=-theta
        while phi < 0:
            phi+=2*m.pi
        while theta < 0:
            theta+=2*m.pi
        self.phi=phi
        self.theta=theta

def readPrescriptionData(zmxPDFile, perturbationFile):
    lines=open(zmxPDFile,'r').readlines()
    length=len(lines)
    i=0 #line number iterator
    while i < length:
        if "GLOBAL VERTEX COORDINATES, ORIENTATIONS, AND ROTATION/OFFSET MATRICES:" in lines[i]:
            perturbationFile=open(perturbationFile,'w')
            k=0 #body command order iterator
            prevSurfNo=-1 #accounts for 4 extra lines
            j=0 #surface number iterator
            while j < len(zmxSurfNos):
                i+=4*(zmxSurfNos[j]-prevSurfNo) #skip directly to next surface
                s=lines[i].split()
                dx=float(s[4])
                rx=float(s[5])
                s=lines[i+1].split()
                dy=float(s[3])
                ry=float(s[4])
                s=lines[i+2].split()
                dz=float(s[3])
                rz=float(s[4])
                #create body command and print to perturbation file
                k=printPerturbation(perturbationFile,BodyCommand(dx,dy,dz,rx,ry,rz,j),k)
                prevSurfNo=zmxSurfNos[j]
                j+=1
            perturbationFile.close() #done
            i=length-1
        i+=1


def printPerturbation(perturbationFile,bodyCommand,i):
    surfID=bodyCommand.surfID
    dx=bodyCommand.dx
    dy=bodyCommand.dy
    dz=bodyCommand.dz
    phi=bodyCommand.phi
    theta=bodyCommand.theta
    if dx != 0:
        perturbationFile.write(names[surfID]+'x '+str(phoSimSurfNos[surfID])+' body 3 random '+str(i)+' 0.0 '+str(dx)+' 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n')
        i+=1
    if dy != 0:
        perturbationFile.write(names[surfID]+'y '+str(phoSimSurfNos[surfID])+' body 4 random '+str(i)+' 0.0 '+str(dy)+' 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n')
        i+=1
    if dz != 0:
        perturbationFile.write(names[surfID]+'z '+str(phoSimSurfNos[surfID])+' body 5 random '+str(i)+' 0.0 '+str(dz)+' 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n')
        i+=1
    if phi != 0:
        perturbationFile.write(names[surfID]+'phi '+str(phoSimSurfNos[surfID])+' body 0 random '+str(i)+' 0.0 '+str(phi)+' 0.0 0.0 0.0 0.0 0.0 '+str(phi)+' 0.0\n')
        i+=1
    if theta != 0:
        perturbationFile.write(names[surfID]+'theta '+str(phoSimSurfNos[surfID])+' body 2 random '+str(i)+' 0.0 '+str(theta)+' 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n')
        i+=1
    return i


def readConfig(configFile):
    #populate list data of surface names and links
    lines=open(configFile,'r').readlines()
    for line in lines:
        s=line.split()
        names.append(s[0])
        zmxSurfNos.append(int(s[1]))
        phoSimSurfNos.append(int(s[2]))


def main():
    configFile=sys.argv[1]
    zmxPDFile=sys.argv[2]
    perturbationFile=sys.argv[3]
    readConfig(configFile)
    readPrescriptionData(zmxPDFile,perturbationFile)
    print('success')


names=[]
zmxSurfNos=[]
phoSimSurfNos=[]

if __name__ == "__main__":
        main()
