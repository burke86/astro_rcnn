#!/usr/bin/env python
##
## @package phosim/tools
## @file phosim_visualizer.py
## @brief python methods called by phosim script
##
## @brief Created by:
## @author Jun Cheng (Purdue), En-Hsin (Purdue)
##
## @brief Modified by:
## @author Colin Burke (Purdue)
##
## @warning This code is not fully validated
## and not ready for full release.  Please
## treat results with caution.
##
## Usage: see README.md file
##
print('------------------------------------------------------------------------------------------')
print('Instrument and Raytrace Visualizer')
print('------------------------------------------------------------------------------------------')

import os,sys
import math
# try to import non-standard modules
try:
    import numpy as np
    from astropy.io import fits
    from mayavi import mlab
except:
    print('Import Error: PhoSim Visualizer failed to import required mayavi python modules (see readme)')


SMALL = 1.0e-12

def convertRotation(ai, aj, ak):
    # convert rotation from phosim (intrinsic zxz) to vtk (extrinsic xyz) 

    # code modified from transformations.py by Christoph Gohlke

    # transformations.py

    # Copyright (c) 2006-2018, Christoph Gohlke
    # Copyright (c) 2006-2018, The Regents of the University of California
    # Produced at the Laboratory for Fluorescence Dynamics
    # All rights reserved.
    #
    # Redistribution and use in source and binary forms, with or without
    # modification, are permitted provided that the following conditions are met:
    #
    # * Redistributions of source code must retain the above copyright
    #   notice, this list of conditions and the following disclaimer.
    # * Redistributions in binary form must reproduce the above copyright
    #   notice, this list of conditions and the following disclaimer in the
    #   documentation and/or other materials provided with the distribution.
    # * Neither the name of the copyright holders nor the names of any
    #   contributors may be used to endorse or promote products derived
    #   from this software without specific prior written permission.
    #
    # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    # AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    # IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    # ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
    # LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    # CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    # SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    # INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    # CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    # ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    # POSSIBILITY OF SUCH DAMAGE.

    axes = 'rzxz'

    # epsilon for testing whether a number is close to zero
    _EPS = np.finfo(float).eps * 4.0

    # axis sequences for Euler angles
    _NEXT_AXIS = [1, 2, 0, 1]

    # map axes strings to/from tuples of inner axis, parity, repetition, frame
    _AXES2TUPLE = {
        'sxyz': (0, 0, 0, 0), 'sxyx': (0, 0, 1, 0), 'sxzy': (0, 1, 0, 0),
        'sxzx': (0, 1, 1, 0), 'syzx': (1, 0, 0, 0), 'syzy': (1, 0, 1, 0),
        'syxz': (1, 1, 0, 0), 'syxy': (1, 1, 1, 0), 'szxy': (2, 0, 0, 0),
        'szxz': (2, 0, 1, 0), 'szyx': (2, 1, 0, 0), 'szyz': (2, 1, 1, 0),
        'rzyx': (0, 0, 0, 1), 'rxyx': (0, 0, 1, 1), 'ryzx': (0, 1, 0, 1),
        'rxzx': (0, 1, 1, 1), 'rxzy': (1, 0, 0, 1), 'ryzy': (1, 0, 1, 1),
        'rzxy': (1, 1, 0, 1), 'ryxy': (1, 1, 1, 1), 'ryxz': (2, 0, 0, 1),
        'rzxz': (2, 0, 1, 1), 'rxyz': (2, 1, 0, 1), 'rzyz': (2, 1, 1, 1)}

    _TUPLE2AXES = dict((v, k) for k, v in _AXES2TUPLE.items())

    # to matrix

    try:
        firstaxis, parity, repetition, frame = _AXES2TUPLE[axes]
    except (AttributeError, KeyError):
        _TUPLE2AXES[axes]  # validation
        firstaxis, parity, repetition, frame = axes

    i = firstaxis
    j = _NEXT_AXIS[i+parity]
    k = _NEXT_AXIS[i-parity+1]

    if frame:
        ai, ak = ak, ai
    if parity:
        ai, aj, ak = -ai, -aj, -ak

    si, sj, sk = math.sin(ai), math.sin(aj), math.sin(ak)
    ci, cj, ck = math.cos(ai), math.cos(aj), math.cos(ak)
    cc, cs = ci*ck, ci*sk
    sc, ss = si*ck, si*sk

    M = np.identity(4)
    if repetition:
        M[i, i] = cj
        M[i, j] = sj*si
        M[i, k] = sj*ci
        M[j, i] = sj*sk
        M[j, j] = -cj*ss+cc
        M[j, k] = -cj*cs-sc
        M[k, i] = -sj*ck
        M[k, j] = cj*sc+cs
        M[k, k] = cj*cc-ss
    else:
        M[i, i] = cj*ck
        M[i, j] = sj*sc-cs
        M[i, k] = sj*cc+ss
        M[j, i] = cj*sk
        M[j, j] = sj*ss+cc
        M[j, k] = sj*cs-sc
        M[k, i] = -sj
        M[k, j] = cj*si
        M[k, k] = cj*ci

    # to euler

    matrix = M

    axes = 'sxyz'
    try:
        firstaxis, parity, repetition, frame = _AXES2TUPLE[axes.lower()]
    except (AttributeError, KeyError):
        _TUPLE2AXES[axes]  # validation
        firstaxis, parity, repetition, frame = axes

    i = firstaxis
    j = _NEXT_AXIS[i+parity]
    k = _NEXT_AXIS[i-parity+1]

    M = np.array(matrix, dtype=np.float64, copy=False)[:3, :3]
    if repetition:
        sy = math.sqrt(M[i, j]*M[i, j] + M[i, k]*M[i, k])
        if sy > _EPS:
            ax = math.atan2( M[i, j],  M[i, k])
            ay = math.atan2( sy,       M[i, i])
            az = math.atan2( M[j, i], -M[k, i])
        else:
            ax = math.atan2(-M[j, k],  M[j, j])
            ay = math.atan2( sy,       M[i, i])
            az = 0.0
    else:
        cy = math.sqrt(M[i, i]*M[i, i] + M[j, i]*M[j, i])
        if cy > _EPS:
            ax = math.atan2( M[k, j],  M[k, k])
            ay = math.atan2(-M[k, i],  cy)
            az = math.atan2( M[j, i],  M[i, i])
        else:
            ax = math.atan2(-M[j, k],  M[j, j])
            ay = math.atan2(-M[k, i],  cy)
            az = 0.0

    if parity:
        ax, ay, az = -ax, -ay, -az
    if frame:
        ax, az = az, ax

    az = 180+az*180/math.pi
    ay = ay*180/math.pi
    ax = ax*180/math.pi
    return ax, ay, az


class photon(object):
    def  __init__(self):
        self.wavelength = 0
        self.color = (0, 0, 0)
        self.xdir = 0
        self.ydir = 0
        self.time = 0
        self.listX = []
        self.listY = []
        self.listZ = []
        self.listLayer = []
        self.opticX = []
        self.opticY = []
        self.opticZ = []

class Control(object):
    def __init__(self, line):
        self.name = line[0]
        self.type = line[3]
        if self.type == '0':
            self.parameter = float(line[13]) # mm (decenter) or rad (tilt)
        else:
            self.parameter = float(line[7]) # mm (decenter) or rad (tilt)
        self.links = [] 
        linkss = line[1].split('|')
        for i in range(len(linkss)):
            self.links.append(int(linkss[i])) 
    def applyControl(self, actor, origin):
        # apply rotation or translation control
        actor.origin = origin
        if self.type == "0":
            alpha, beta, gamma = convertRotation(self.parameter, 0, 0)
            actor.orientation = [actor.orientation[0], actor.orientation[1], actor.orientation[2]+gamma]
        elif self.type == "1":
            alpha, beta, gamma = convertRotation(0, 0, self.parameter)
            actor.orientation = [actor.orientation[0], actor.orientation[1]+beta, actor.orientation[2]]
        elif self.type == "2":
            alpha, beta, gamma = convertRotation(0, -self.parameter, 0)
            actor.orientation = [actor.orientation[0]+alpha, actor.orientation[1], actor.orientation[2]]
        elif self.type == "3":
            actor.position = [actor.position[0]+self.parameter, actor.position[1], actor.position[2]]
        elif self.type == "4":
            actor.position = [actor.position[0], actor.position[1]+self.parameter, actor.position[2]]
        elif self.type == "5":
            actor.position = [actor.position[0], actor.position[1], actor.position[2]+self.parameter]
        
class Surface(object): 
    # surface constructor: read one line for a specific surface
    def __init__(self, line, surfaceNumber): 
        self.name = line[0]
        self.type = line[1]
        self.curvature = float(line[2]) # mm
        self.thickness = float(line[3]) # mm
        self.outerRadius = float(line[4]) # mm
        self.innerRadius = float(line[5]) # mm
        self.conic = float(line[6]) # no units
        self.aspheric3 = float(line[7])*1e3 # mm
        self.aspheric4 = float(line[8])*1e3 # mm
        self.aspheric5 = float(line[9])*1e3 # mm
        self.aspheric6 = float(line[10])*1e3 # mm
        self.aspheric7 = float(line[11])*1e3 # mm
        self.aspheric8 = float(line[12])*1e3 # mm
        self.aspheric9 = float(line[13])*1e3 # mm
        self.aspheric10 = float(line[14])*1e3 # mm
        self.coatingFile = line[15]
        self.mediaFile = line[16]
        self.surfaceNumber = surfaceNumber
        self.actor = None
    def plotSurface(self, globalZ, controls):
        # plot optical surfaces
        if(self.type!="det"):
            origin = [0, 0, globalZ]
            r, theta = np.mgrid[self.innerRadius:self.outerRadius:100j, -np.pi:np.pi:100j]
            # sag equation
            if abs(self.curvature) < SMALL:
                z = globalZ + (r**3*self.aspheric3 + r**4*self.aspheric4 + r**5*self.aspheric5 + r**6*self.aspheric6 
                    + r**7*self.aspheric7 + r**8*self.aspheric8 + r**9*self.aspheric9 + r**10*self.aspheric10)
            else:
                z = globalZ + (r**2/(self.curvature*(1 + np.sqrt(1 - (1 + self.conic)*r**2/(self.curvature**2))))
                    + r**3*self.aspheric3 + r**4*self.aspheric4 + r**5*self.aspheric5 + r**6*self.aspheric6 
                    + r**7*self.aspheric7 + r**8*self.aspheric8 + r**9*self.aspheric9 + r**10*self.aspheric10)
            x = r*np.cos(theta)
            y = r*np.sin(theta)
            r = np.sqrt(x**2 + y**2)
            if self.type == "filter": 
                mesh = mlab.mesh(x,y,z, opacity=0.5, color = (1,1,0))
            elif self.type == "none":
                mesh = mlab.mesh(x,y,z, opacity=0.05, color=(0.2,0.7,0.9))
            elif self.type == "lens":
                mesh = mlab.mesh(x,y,z, opacity=0.5, color=(0.2,0.7,0.9))
            else:
                mesh = mlab.mesh(x,y,z, opacity=1.0, color=(0.9,0.9,0.9), colormap="Pastel2")
            
            for control in controls:
                for link in control.links:
                    if link == self.surfaceNumber:
                        control.applyControl(mesh.actor.actor, origin)
            self.actor = mesh.actor.actor
                
class Chip(object): 
        def __init__(self,line, surfaceNumber): 
            self.name = line[0]
            self.centerX = float(line[1])   # microns
            self.centerY = float(line[2])   # microns
            self.pixelSize = float(line[3]) # microns
            self.numX = float(line[4]) # pixels
            self.numY = float(line[5]) # pixels
            self.thickness = float(line[10]) # microns
            self.phi = float(line[12])*math.pi/180 # degrees
            self.theta = -float(line[13])*math.pi/180 # degrees
            self.psi = float(line[14])*math.pi/180 # degrees

            self.halfX = self.pixelSize *self.numX/2
            self.halfY = self.pixelSize *self.numY/2
            self.halfZ = self.thickness/2
            self.surfaceNumber = surfaceNumber #surface number (used to link control.txt)

        def plotChip(self, pos, controls):
            # plot chips in focal plane
            x, y = np.mgrid[-self.halfX:self.halfX:2j, -self.halfY:self.halfY:2j]
            x = (x+self.centerX)/1000
            y = (y+self.centerY)/1000
            z = x-x+pos
            mesh = mlab.mesh(x,y,z, opacity = 0.95, color=(0.1,0.1,0.1))
            for control in controls:
                for link in control.links:
                    if link == self.surfaceNumber:
                        control.applyControl(mesh.actor.actor, [0, 0, pos])
            # additional focalplanelayout.txt rotations
            if self.phi != 0 or self.theta != 0 or self.psi != 0:
                mesh.actor.actor.origin = [self.centerX/1000, self.centerY/1000, pos]
                alpha, beta, gamma = convertRotation(self.phi, self.theta, self.psi)
                mesh.actor.actor.orientation = [mesh.actor.actor.orientation[0]+alpha, mesh.actor.actor.orientation[1]+beta, mesh.actor.actor.orientation[2]+gamma]

def readArgs(instrDir, eventFiles):
    # choose first eventFile (they must have the same characteristics) and get file names
    eventFile = eventFiles.split()[0]
    instrument = instrDir.split('/')[-2]
    s = eventFile.split('/')[-1].replace(instrument,'').split('_')[1:]
    if s[0] != 'r':
        print('Not an eventfile!')
        sys.exit()
    observationID = s[1]
    filt = s[2][1:]

    opticsFile = instrDir + 'optics_'+filt+'.txt'
    focalplaneFile = instrDir + 'focalplanelayout.txt'
    segmentationFile = instrDir + 'segmentation.txt'
    perturbationFile = instrDir + 'perturbation.txt'
    spiderFile = instrDir + 'spider.txt'

    return opticsFile, focalplaneFile, segmentationFile, perturbationFile, spiderFile

def readMultpleEvents(eventFiles, maxZ, minZ):
    # read list of eventfiles
    for eventFile in eventFiles.split():
        readEvents(eventFile, maxZ, minZ)


def readEvents(eventFits, maxZ, minZ):
    # read events in eventfile by layer
    if os.path.exists(eventFits):

        MAXLINE = 1000
        print("Plotting photon rays (first "+str(MAXLINE)+" lines).")
        f = fits.open(eventFits)
        event=f[1].data
        f.close()
        xpos = []
        ypos = []
        zpos = []
        wavelength = []
        numLine = len(event.field(0))
        if numLine > MAXLINE:
            numLine = MAXLINE
        numPhoton = 0
        photonList = []
        pre_p = photon()
        for i in range(0, numLine):
            layer = event.field(3)[i]
            if layer == 0: # photon properties
                if i != 0:
                    photonList.append(pre_p)
                numPhoton += 1
                p = photon()
                p.wavelength = event.field(2)[i]
                p.color = photonColor(p.wavelength)
                p.xdir = event.field(0)[i]
                p.ydir = event.field(1)[i]
            elif layer == 1: # time
                p.time = event.field(0)[i]
            elif layer == 99: # initial point
                z = maxZ - minZ
                dz = event.field(2)[i] - z
                p.listX.append(event.field(0)[i] + dz*p.xdir)
                p.listY.append(event.field(1)[i] + dz*p.ydir)
                p.listZ.append(z)
                p.listLayer.append(layer)
            elif layer > 199 and layer < 300: # telescope
                p.listX.append(event.field(0)[i])
                p.listY.append(event.field(1)[i])
                p.listZ.append(event.field(2)[i])
                p.listLayer.append(layer)

            pre_p = p

        photonList.append(p)

        for i in range(len(photonList)):
            photoni = photonList[i]
            for n in range(len(photoni.listZ)):
                mlab.plot3d(photoni.listX[n:2+n], photoni.listY[n:2+n], photoni.listZ[n:2+n], color = photoni.color, opacity = 0.2, tube_radius = None)

def readBodyCommands(perturbationFile):
    # read perturbation.txt and return list of body commands
    print('Plotting optics.')
    controls = []
    if os.path.exists(perturbationFile):
        for line in open(perturbationFile).readlines():
            if line[0] != "#" and line.strip() != "":
                controls.append(Control(line.split()))
    return controls

def readOptics(opticsFile, controls):
    # read optics.txt and plot
    surfaces = []
    surfno = 0
    globalZ = 0.0
    globalZControl = 0.0
    maxZ = -100000.0
    # gather surfaces from optics file
    for line in open(opticsFile).readlines():
        if line[0] != "#":
            surfaces.append(Surface(line.split(), surfno))
            surfno+=1
    # loop through surfaces and plot
    for i in range(len(surfaces)):
        globalZ += surfaces[i].thickness
        # check for where to start rays
        if globalZ > maxZ:
            maxZ = globalZ
        surfaces[i].plotSurface(globalZ, controls)
        # double check maxZ after controls are applied
        if surfaces[0].curvature < 0: # refracting telescope
            maxZ = surfaces[0].curvature*2 
        elif surfaces[i].actor != None:
            globalZControl += surfaces[i].actor.position[2]
            if globalZControl > maxZ:
                maxZ = globalZControl
    minZ = surfaces[0].actor.position[2]
    return maxZ, minZ, globalZ, surfno - 1

def readChips(focalplaneFile, detPosition, surfno, controls):
    # read focalplanelayout.txt and plot
    print("Plotting focal plane.")
    chips = []
    for line in open(focalplaneFile).readlines():
        if line[0] != "#":
            chips.append(Chip(line.split(), surfno))
    for i in range(len(chips)):
        chips[i].plotChip(detPosition, controls)

def photonColor(wavelength):
    # approximate wavelength to rgb color correction
    # wavelength in nm
    # follows www.physics.sfasu.edu/astro/color/spectra.html (Dan Bruton)

    wavelength *= 1000 # um
    gamma = 0.8 # intensity parameter

    if wavelength < 380: # UV
        r = 0.5
        g = 0.0
        b = 0.5
    if wavelength >= 380 and wavelength <= 440: # violet
        r = ((-(wavelength - 440)/(440 - 380))*(0.3 + 0.7*(wavelength - 380)/(440 - 380)))**gamma
        g = 0.0
        b = (1.0*0.3 + 0.7*(wavelength - 380)/(440 - 380))**gamma
    elif wavelength >= 440 and wavelength <= 490: # blue
        r = 0.0
        g = ((wavelength - 440)/(490 - 440))**gamma
        b = 1.0
    elif wavelength >= 490 and wavelength <= 510: # green
        r = 0.0
        g = 1.0
        b = (-(wavelength - 510)/(510 - 490))**gamma
    elif wavelength >= 510 and wavelength <= 580: # yellow
        r = ((wavelength - 510)/(580 - 510))**gamma
        g = 1.0
        b = 0.0
    elif wavelength >= 580 and wavelength <= 645: # orange
        r = 1.0
        g = (-(wavelength - 645)/(645 - 580))**gamma
        b = 0.0
    elif wavelength >= 645 and wavelength <= 750: # red
        r = (0.3 + 0.7*(750 - wavelength)/(750 - 645))**gamma
        g = 0.0
        b = 0.0
    else: # IR
        r = 1.0
        g = 0.0
        b = 0.0

    return (r, g, b)

def main():
    # read arguments
    if len(sys.argv) < 3:
        print('Not enough arguments! Requires: <ISC dir> <list of eventfile paths>')
        sys.exit()
    instrDir = sys.argv[1]
    eventFiles = sys.argv[2]
    if instrDir[-1] != '/':
        instrDir += '/'
    opticsFile, focalplaneFile, segmentationFile, perturbationFile, spiderFile = readArgs(instrDir, eventFiles)

    # set figure paramaters
    fig = mlab.figure(bgcolor = (0.99, 0.99, 0.99), size = (512, 512))

    # read data and plot
    controls = readBodyCommands(perturbationFile)
    maxZ, minZ, detPosition, surfno = readOptics(opticsFile, controls)
    readChips(focalplaneFile, detPosition, surfno, controls)
    readMultpleEvents(eventFiles, maxZ, minZ)

    mlab.show()
    
if __name__ == "__main__":
        main()
