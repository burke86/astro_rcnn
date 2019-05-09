#!/usr/bin/env python
import subprocess
import shutil, sys, os, glob, optparse
import re

class Validation:
    def __init__(self,opt):
        self.outputDir=opt.outputDir
        self.grid=opt.grid
        self.validationScript='validation_pipeline'
        self.tlogFile=opt.tlog
        self.flat=opt.flat
        self.task=opt.task

        if self.task=='all':
            self.unit=True
            self.speed=True
        elif self.task=='none':
            self.unit=False
            self.speed=False
        else:
            lstr = self.task.split('|')
            if 'unit' in lstr:
                self.unit=True
            else:
                self.unit=False
            if 'speed' in lstr:
                self.speed=True
            else:
                self.speed=False


        if self.grid=='PBS':
            runProgram("source /etc/profile.d/modules.sh; module load gcc/4.7.2; module load python/2.7; python -V; gcc -v")
            runProgram("git pull; make")
        if self.outputDir=='new':
            version=open('bin/version').readlines()
            self.outputDir='validation_'+version[1].split()[0]+'/'
        if not os.path.exists(self.outputDir):
            os.makedirs(self.outputDir)

        self.tlogID=[]
        self.tlogHH=[]
        self.tlogMM=[]

        if self.tlogFile is not None:
            for line in open(self.tlogFile):
                lstr=line.split()
                self.tlogID.append(lstr[0])
                estimateTime=float(lstr[1])
                actualTime=float(lstr[2])
                if estimateTime-actualTime<1.0 or estimateTime/actualTime<2:
                    self.tlogHH.append(int(estimateTime))
                    self.tlogMM.append(int((estimateTime-int(estimateTime))*60))
                else:
                    estimateTime=actualTime*2
                    self.tlogHH.append(int(estimateTime))
                    self.tlogMM.append(int((estimateTime-int(estimateTime))*60))

    def setup(self):  #read validation_pipeline script
        self.jobID=[]
        taskID=[]
        subTaskNumber=[]
        self.catalog=[]
        self.command=[]
        self.args=[]
        self.runLocally=[]
        self.localCommand=[]
        self.runFlag=[]

        for line in open(self.validationScript):
            if 'phosim' in line and 'catalog' in line and '4A' not in line:
                s=re.split(';|\(|\)|>|\n',line)
                for i in range(len(s)):
                    if 'phosim' in s[i]:
                        ss=s[i].split()
                        cat=ss[1]
                        taskid=cat.split('_')[1]
                        self.catalog.append(cat)
                        idx=[k for k,x in enumerate(taskID) if x==taskid]
                        if len(idx)==0:
                            taskID.append(taskid)
                            subTaskNumber.append(0)
                            jid='%s_%02d' % (taskid,0)
                            self.jobID.append(jid)
                        else:
                            subTaskNumber[idx[0]]+=1
                            jid='%s_%02d' % (taskid,subTaskNumber[idx[0]])
                            self.jobID.append(jid)
                        args=''
                        j=2
                        while j<len(ss):
                            if ss[j]=='-c':
                                self.command.append(ss[j+1])
                                j+=1
                            elif ss[j]=='-o' or ss[j]=='-p': #ignore outputDir and Nprocessor setting in validation_pipeline
                                j+=1
                            else:
                                args=args+' '+ss[j]
                            j+=1
                        self.args.append(args)
                        self.runLocally.append(False)
                        self.localCommand.append(None)
                        if self.task=='all':
                            self.runFlag.append(True)
                        elif self.task=='none':
                            self.runFlag.append(False)
                        else:
                            self.runFlag.append(False)
                            lstr = self.task.split('|')
                            for j in range(len(lstr)):
                                if lstr[j]==taskid:
                                    self.runFlag[-1]=True
                                    break
                        break
            if 'mv' in line:  #jobs have to be run locally
                self.runLocally[-1]=True
                self.localCommand[-1]=line.replace('validation/',self.outputDir).split('\n')[0]

        #for i in range(len(self.catalog)):
        #    print self.jobID[i], self.catalog[i], self.command[j]

    def findCPUtime(self,jobID):
        idx=[k for k,x in enumerate(self.tlogID) if x==jobID]
        if len(idx)==0:
            hh=1
            mm=0
        else:
            hh=self.tlogHH[idx[0]]
            mm=self.tlogMM[idx[0]]
        return hh, mm

    def submitPBSjob(self,i):
        hh, mm=self.findCPUtime(self.jobID[i])
        sub=self.jobID[i]+'.sub'
        submitfile=open(sub,'w')
        submitfile.write('#!/bin/sh -l\n')
        submitfile.write('#PBS -l nodes=1:ppn=1\n')
        submitfile.write('#PBS -l walltime=%02d:%02d:00\n' % (hh,mm))
        submitfile.write('module load gcc/4.7.2\n')
        submitfile.write('module load python/2.7\n')
        submitfile.write('cd $PBS_O_WORKDIR\n')
        submitfile.write('unset DISPLAY\n')
        submitfile.write('python phosim.py '+self.catalog[i]+' -c '+self.command[i]+' -o '+self.outputDir+' '+self.args[i]+'\n')
        submitfile.close()
        if hh>=4:
            command='qsub -q peters11 '+sub
        else:
            command='qsub '+sub
        #print command
        runProgram(command)

    def runUnitTest(self):
        outputDir=self.outputDir
        runProgram('validation/unittest > '+outputDir+'unittest_output')
        os.chdir('source/ancillary/')
        runProgram('./test_phosim_parser < test.pars >> ../../'+outputDir+'unittest_output')
        runProgram('./test_readtext < test.txt >> ../../'+outputDir+'unittest_output')
        os.chdir('../../')

    def runSpeedTest(self):
        outputDir=self.outputDir
        if self.grid=='PBS':
            hh, mm=self.findCPUtime('4A')
            submitfile=open('4A.sub','w')
            submitfile.write('#!/bin/sh -l\n')
            submitfile.write('#PBS -l nodes=1:ppn=1\n')
            submitfile.write('#PBS -l walltime=%02d:%02d:00\n' % (hh,mm))
            submitfile.write('module load gcc/4.7.2\n')
            submitfile.write('module load python/2.7\n')
            submitfile.write('cd $PBS_O_WORKDIR\n')
            submitfile.write('unset DISPLAY\n')
            submitfile.write('date > '+outputDir+'speed_out_0\n')
            submitfile.write('python phosim.py validation/validation_4A_00_catalog -c validation/validation_4A_00_commands -o '+outputDir+'\n')
            submitfile.write('date > '+outputDir+'speed_out_1\n')
            submitfile.write('python phosim.py validation/validation_4A_01_catalog -c validation/validation_4A_01_commands -o '+outputDir+' >& '+outputDir+'speed_out\n')
            submitfile.write('date > '+outputDir+'speed_out_2\n')
            submitfile.close()
            runProgram('qsub -q peters11 4A.sub')
        else:
            runProgram('date > '+outputDir+'speed_out_0')
            runProgram('python phosim.py validation/validation_4A_00_catalog -c validation/validation_4A_00_commands -o '+outputDir)
            runProgram('date > '+outputDir+'speed_out_1')
            runProgram('python phosim.py validation/validation_4A_01_catalog -c validation/validation_4A_01_commands -o '+outputDir+' >& '+outputDir+'speed_out')
            runProgram('date > '+outputDir+'speed_out_2')

    def generateFlats(self):
        outputDir=self.outputDir
        filters=['u','g','r','i','z','y']
        if self.grid=='PBS':
            for f in range(6):
                #hh, mm=self.findCPUtime('F_'+str(f))
                hh=144
                mm=0
                sub='F_'+str(f)+'.sub'
                submitfile=open(sub,'w')
                submitfile.write('#!/bin/sh -l\n')
                submitfile.write('#PBS -l nodes=1:ppn=1\n')
                submitfile.write('#PBS -l walltime=%02d:%02d:00\n' % (hh,mm))
                submitfile.write('module load gcc/4.7.2\n')
                submitfile.write('module load python/2.7\n')
                submitfile.write('cd $PBS_O_WORKDIR\n')
                submitfile.write('unset DISPLAY\n')
                submitfile.write('time python phosim.py examples/flats/flat'+filters[f]+'_instcat_1 -s R22_S11 -o '+outputDir+'\n')
                submitfile.close()
                runProgram('qsub -q peters11 '+sub)
        else:
            for f in range(6):
                runProgram('time python phosim.py examples/flats/flat'+filters[f]+'_instcat_0 -s R22_S11 -o '+outputDir)


    def run(self):
        for i in range(len(self.catalog)):
            if self.runFlag[i]:
                if self.grid=='PBS' and not self.runLocally[i]:
                    self.submitPBSjob(i)
                else:
                    runProgram('python phosim.py '+self.catalog[i]+' -c '+self.command[i]+' -o '+self.outputDir+' '+self.args[i])
                    if self.localCommand[i] is not None:
                        runProgram(self.localCommand[i])
        if self.unit:
            self.runUnitTest()
        if self.speed:
            self.runSpeedTest()
        if self.flat==1:
            self.generateFlats()


def runProgram(command):
    if subprocess.call(command, shell=True) != 0:
        raise RuntimeError("Error running %s" % command)

def removeFile(filename):
    try:
         os.remove(filename)
    except OSError:
         pass




def main():
    parser = optparse.OptionParser()
    parser.add_option('-o','--output',dest="outputDir",default='validation/')
    parser.add_option('-T','--tlog',dest="tlog",default=None)
    parser.add_option('-g','--grid',dest="grid",default="no")
    parser.add_option('-f','--flat',dest="flat",default=0,type="int")
    parser.add_option('-t','--task',dest="task",default="all")
    opt, remainder = parser.parse_args(sys.argv[1:])

    validation=Validation(opt)
    validation.setup()
    validation.run()


if __name__ == "__main__":
    main()

