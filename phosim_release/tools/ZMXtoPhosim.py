#!/usr/bin/env python

"""
  @package phosim
  @file ZMXtoPhosim.py
  @brief python script to convert a ZEMAX file to phosim optics_*.txt
 
  @brief Created by:
  @author En-Hsin Peng (Purdue)
 
  @warning This code is not fully validated
  and not ready for full release.  Please
  treat results with caution.

Usage: python ZMXtoPhosim.py yourZemax.ZMX input.lis

       input.lis is a user provided file which contains
       Column 1: Zemax surface number 
       Column 2: Surface name (M1, L1, etc) 
       Column 3: Surface type (mirror, lens, filter, det, grating) 
       Column 4: Coating file
       Column 5: Medium file

       e.g.,
       23   M1    mirror  m1_protAl_Ideal.txt  air
       26   M2    mirror  m2_protAl_Ideal.txt  air
       29   M3    mirror  m3_protAl_Ideal.txt  air
       30   none  none    none                 air
       31   L1    lens    lenses.txt           silica_dispersion.txt
       32   L1E   lens    lenses.txt           air
       33   L2    lens    lenses.txt           silica_dispersion.txt
       34   L2E   lens    lenses.txt           air
       35   F     filter  filter_x.txt         silica_dispersion.txt
       36   FE    filter  none                 air
       37   L3    lens    lenses.txt           silica_dispersion.txt
       38   L3E   lens    lenses.txt           air
       39   D     det     detectorar.txt       air

Notes: 1. This script assumes ZEMAX in mm. 
       2. It will overwrite existing optics_*.txt in the current folder. 

"""

import sys, os, commands

class Surface(object):
    def __init__(self,curv,disz,z,rout,rin,coni,an):
        self.curv=curv
        self.disz=disz
        self.z=z
        self.rout=rout
        self.rin=rin
        self.coni=coni
        self.an=an

def findSurface(zmx,surf,surf0,flt,zprev=0.0):
    an=[0.0 for i in range(16)]
    coni=0.0
    curv=0.0
    disz=[0.0 for i in range(100)]
    pzup=[[-1.0,-1.0,-1.0] for i in range(100)]
    rout=0.0
    rin=0.0
    readCurv=False
    #read DISZ, update THIC, CRVT
    for line in open(zmx).readlines():
        if 'SURF ' in line:
            s=int(float(line.split()[1]))
        elif 'DISZ' in line:
            if line.split()[1]!='INFINITY':
                disz[s]=float(line.split()[1])
        elif 'THIC' in line:
            l=line.split()
            if float(l[2])-1==flt:
                disz[int(float(l[1]))]=float(l[3])
        elif 'CRVT' in line:
            l=line.split()
            if float(l[2])-1==flt and float(l[1])==surf:
                curv=-1/float(l[3])
                readCurv=True
        elif 'PZUP' in line:
            l=line.split()
            pzup[s][0]=float(l[1])
            pzup[s][1]=float(l[2])
            pzup[s][2]=float(l[3])

    #PZUP
    for i in range(surf):
        if pzup[i][0]>0:
            disz[i]=disz[int(pzup[i][0])]*pzup[i][1]+pzup[i][2]


    z=0.0
    if surf>surf0:
        for i in range(surf):
            if i==surf0:
                z0=z
            z-=disz[i]
        z-=z0

    #read asphere
    found=False
    for line in open(zmx).readlines():
        if 'SURF '+str(surf) in line:
            found=True
            continue
        if found:
            if 'CURV' in line and readCurv==False:
                if float(line.split()[1])!=0:
                    curv=-1/float(line.split()[1])
            elif 'CLAP' in line or 'FLAP' in line:
                rout=float(line.split()[2])
                rin=float(line.split()[1])
            elif 'DIAM' in line:
                rout=float(line.split()[1])
                rin=0.0
            elif 'CONI' in line:
                coni=float(line.split()[1])
            elif 'PARM' in line:  #parm 1: a2, parm 2: a4
                lstr=line.split()
                an[int(float(lstr[1]))*2-1]=float(lstr[2])/1e3
            elif 'SURF' in line or 'CONF' in line:
                    break

    surface=Surface(curv,z-zprev,z,rout,rin,coni,an)
    return surface


def printOptics(output,surface,name,typ,coating,medium,flt):
    out=open(output,'a')
    if typ=='det' and surface.rout==0:
        surface.rout=400.0
    out.write('%-7s %7s %8.1f %10.4f %10.4f %6.1f %6.3f '  % (name,typ,surface.curv,surface.disz,surface.rout,surface.rin,surface.coni))
    for i in range(8):
        out.write('%10.3e ' % (surface.an[i+2]))
    if coating!='none' and typ=='filter':
        coating='filter_'+str(flt)+'.txt'
    out.write('%s %s\n' % (coating,medium))
    out.close()


zmxFile=sys.argv[1]
inputList=sys.argv[2]
fltNum=int(float(commands.getoutput('grep MNUM '+zmxFile+' | awk \'{print $2}\'')))

zmxSurface=[]
name=[]
typ=[]
coating=[]
medium=[]
for line in open(inputList).readlines():
    l=line.split()
    zmxSurface.append(int(float(l[0])))
    name.append(l[1])
    typ.append(l[2])
    coating.append(l[3])
    medium.append(l[4])


for flt in range(fltNum):
    output='optics_'+str(flt)+'.txt'
    try:
        os.remove(output)
    except OSError:
        pass

    zprev=0.0
    for i in range(len(name)):
        surface=findSurface(zmxFile,zmxSurface[i],zmxSurface[0],flt,zprev)
        zprev=surface.z
        printOptics(output,surface,name[i],typ[i],coating[i],medium[i],flt)
