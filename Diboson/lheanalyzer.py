#!/usr/bin/env python
import gzip
class LHEFileFormatError(Exception):
    def __init__(self, filename,line):
        self.line = line
        self.filename = filename
    def __str__(self):
        return repr("Error in LHE File "+self.filename+" in line "+str(self.line))

class Particle:
    """Describes a single particle"""
    def __init__(self,initstr):
        """Constuctor call with line from lhe file"""
        ls=initstr.split()
        self.pdgId,self.status,self.mother, self.color,self.momentum,self.lifetime,self.spin=int(ls[0]),int(ls[1]),[int(ls[2]),int(ls[3])],[int(ls[4]),int(ls[5])],[float(ls[6]),float(ls[7]),float(ls[8]),float(ls[9]),float(ls[10])],float(ls[11]),float(ls[12])
    def getTLorentzVector(self):
        """Returns the ROOT TLorentzVector of the particle"""
        import ROOT
        return ROOT.TLorentzVector(self.momentum[0],self.momentum[1],self.momentum[2],self.momentum[3])
    px = property(lambda self: self.momentum[0])
    py = property(lambda self: self.momentum[1])
    pz = property(lambda self: self.momentum[2])
    energy = property(lambda self: self.momentum[3])
    mass = property(lambda self: self.momentum[4])
    pt = property(lambda self: (self.px**2+self.py**2)**0.5)
    tLorentzVector = property(getTLorentzVector)
    eta = property(lambda self: self.tLorentzVector.Eta())
    phi = property(lambda self: self.tLorentzVector.Phi())
       
class Event:
    def __init__(self,initstr):
        """Constuctor call with line from lhe file"""
        ls=initstr.split()
        self.nParticles,self.processId,self.weight,self.scale,self.QEDCoupling,self.QCDCoupling=int(ls[0]),int(ls[1]),float(ls[2]),float(ls[3]),float(ls[4]),float(ls[5])
        self.particles=[]
        self.rwgts={}
    def addParticle(self,particle):
        """adds a particle to the event"""
        self.particles.append(particle)

class Process:
    def __init__(self,initstr):
        ls=initstr.split()
        self.crossSection,self.crossSectionUncertainty,self.maxWeight,self.id=float(ls[0]),float(ls[1]),float(ls[2]),int(ls[3])
        

class LHEFile:
    def __init__(self, filename):
        self.filename=filename
        if filename.endswith('.gz'):
            self.fp = gzip.open(filename, 'r')
        else:           
            self.fp = open(filename,"r")
        self.lineCounter=0
        for line in self.fp:
            self.lineCounter+=1
            if line.strip()=="<init>": break
    def __iter__(self):
        return self
    def next(self):
        for line in self.fp:
            self.lineCounter+=1
            if line[0]!="#":
                return line
                break
        self.fp.close()
                
class LHEAnalysis:
    """Iterator for looping over a sequence backwards."""
    def __init__(self, filename):
        self.lhefile=LHEFile(filename)
	line=self.lhefile.next()
	ls = line.split()
	try:
            self.beamId ,self.beamEnergy,self.PDFAuthor,self.PDFSet,self.weightSwitch,self.nProcesses=[int(ls[0]),int(ls[1])],[float(ls[2]),float(ls[3])],[int(ls[4]),int(ls[5])],[int(ls[6]),int(ls[7])],int(ls[8]),int(ls[9])
        except (ValueError, IndexError):
            print ls
            raise LHEFileFormatError(self.lhefile.filename,self.lhefile.lineCounter)
	self.processes=[]
	for i in range(self.nProcesses):
            line=self.lhefile.next()
            try:
                self.processes.append(Process(line))
            except (ValueError, IndexError):
                print ls
                raise LHEFileFormatError(self.lhefile.filename,self.lhefile.lineCounter)
      
    def __iter__(self):
        return self
    def next(self):
        beginevent=False
        endevent=False
        for line in self.lhefile:
            if line==None:
                raise StopIteration
                return
            if line.strip()=="<event>": 
                beginevent=True
                break
        if beginevent is False:
            raise StopIteration
        initline=self.lhefile.next()
        event=Event(initline)
        for line in self.lhefile:
            if line.strip()=="</event>":
                endevent=True
                break
            if line.strip()[0]=="<":
                break
            try:
                particle=Particle(line)
            except (ValueError, IndexError):
                print line
                raise LHEFileFormatError(self.lhefile.filename,self.lhefile.lineCounter)
            event.addParticle(particle)
        if not endevent:
            beginweights=False
            for line in self.lhefile:
                if line.strip()=="</event>":
                    break
                if line.strip()=="<rwgt>":
                    beginweights=True
                    break
            if beginweights:
                for line in self.lhefile:
                    if line.strip()=="</rwgt>" or line.strip()=="</event>":
                        break
                    else:
                        ls=line.split()
                        assert ls[0]=='<wgt'
                        assert ls[3]=='</wgt>'
                        assert ls[1].split('=')[0]=='id'
                        rwgtid=ls[1].split('=')[1][1:-2]
                        rwgtval=float(ls[2])
                        event.rwgts[rwgtid]=rwgtval
                        rwgtid=line[line.find('<')+1:line.find('>')].split('=')[1][1:-1]
                        rwgtval=float(line[line.find('>')+1:line.find('<',line.find('<')+1)])
                        event.rwgts[rwgtid]=rwgtval

        return event
    totalCrossSection = property(lambda self: sum(p.crossSection for p in self.processes))
    
class LorentzVector:
    def __init__(self,momentum):
        self.px=momentum[0]
        self.py=momentum[1]
        self.pz=momentum[2]
        self.energy=momentum[3]
    def add(self,p):
        self.energy+=p.energy
        self.px+=p.px
        self.py+=p.py
        self.pz+=p.pz
    def invariantMass(self):
        return (self.energy**2-self.px**2-self.py**2-self.pz**2)**0.5
    def momentum(self):
        return ( [px,py,pz] )
        
    
def invariantMass(p0,*particles):
    l0=LorentzVector(p0.momentum)
    for p in particles:
        l0.add(LorentzVector(p.momentum))
    return l0.invariantMass()
def transverseMass(p0,p1):
    import ROOT
    from math import cos
    l0=p0.tLorentzVector
    l1=p1.tLorentzVector
    return (l0.Et()**2+l1.Et()**2-2*l0.Et()*l1.Et()*cos(l0.DeltaPhi(l1)))**0.5


