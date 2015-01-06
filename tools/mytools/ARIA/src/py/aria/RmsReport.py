"""
Authors: Bardiaux Benjamin
         Structural Bioinformatics, Institut Pasteur, Paris
        
         Copyright (C) Benjamin Bardiaux
         No warranty implied or expressed.
         All rights reserved.

$Author: bardiaux $
$Revision: 1.1.1.1 $
$Date: 2010/03/23 15:27:24 $
"""

import os, sys

from aria.ariabase import *
from numpy import *

RMS_TEXT_REPORT = 'noe_restraints.rms'
RMS_PS_REPORT = 'rms_analysis.ps'

RMS_PS_MAP_LEGEND = 'rms_legend.eps'
RMS_PS_MAP = 'rms_2D_residue_map.eps'
RMS_PS_PROFILE = 'rms_residue_profile.eps'

_use_matplotlib = 0 
try:
    import matplotlib
    _use_matplotlib = 1
    matplotlib.use('PS', warn=False)
except:
    pass
    
class RmsReport(AriaBaseClass):



    def __init__(self, peaks, iteration_n, text_path, graphics_path):

        AriaBaseClass.__init__(self)
        self.peaks = peaks
        self.graphics_path = graphics_path
        self.text_report = os.path.join(text_path, RMS_TEXT_REPORT)
        self.ps_report =  os.path.join(graphics_path, RMS_PS_REPORT)
        self.iteration_n = iteration_n

        self.dimer = None
        
    def getOverallRms(self, peaks):

        if not len(peaks):
            return 0.

        viol = power(array([ p.analysis.getUpperBoundViolation()[0] for p in peaks]), 2)
        rms = sum(viol) / len(viol)
        rms = sqrt(rms)        

        return rms

    def getRmsPerResidue(self, peaks):
 
        violations = {}
        weights = {}

        for p in peaks:

            target_dist = p.getUpperBound()
            for c in p.getActiveContributions():
                r =  c.getSpinSystems()[0].getAtoms()[0].getResidue()
                r1 = r.getNumber()
##                 s1 = r.getChain().getSegid()

                r = c.getSpinSystems()[1].getAtoms()[0].getResidue()
                r2 = r.getNumber()
                
##                 s2 = r.getChain().getSegid()
                
                d = c.getAverageDistance()[0]
                w = c.getWeight()
                viol = d - target_dist

    ##             if w >0:
    ##                 print p.getId(), p.analysis.getUpperBoundViolation()[0], d, target_dist, (viol), w
                
                violations.setdefault(r1, [])
                weights.setdefault(r1, [])
                violations[r1].append(viol)
                weights[r1].append(w)
                

                
                if r1 <> r2:
                    

                    
                    violations.setdefault(r2, [])
                    weights.setdefault(r2, [])
                    violations[r2].append(viol)
                    weights[r2].append(w)

        rms_per_resid = {}

        for r,v in violations.items():

            viol = array(v) 
            viol = power(viol * greater(viol, 0), 2) * array(weights[r])
            rms  = sum(viol) / sum(weights[r])
            rms = sqrt(rms)
            rms_per_resid[r] = rms
                


        
        return rms_per_resid

    def getRmsPerResiduePair(self, peaks):


        residues = {}
        weights = {}
        
        for p in peaks:

            target_dist = p.getUpperBound()
            for c in p.getActiveContributions():
                r1 = c.getSpinSystems()[0].getAtoms()[0].getResidue()
                #r1 = r.getNumber()
                #s1 = r.getChain().getSegid()

                r2 = c.getSpinSystems()[1].getAtoms()[0].getResidue()
                #r2 = r.getNumber()
                #s2 = r.getChain().getSegid()

                # check if we have a dimer
                if not self.dimer:
                    if r1.getChain().getSegid() <> r2.getChain().getSegid():
                        self.dimer = 1
                    else:
                        self.dimer = 0

                d = c.getAverageDistance()[0]
                w = c.getWeight()
                viol = d - target_dist

                    
                residues.setdefault(r1, {})
                weights.setdefault(r1, {})
                residues[r1].setdefault(r2, [])
                weights[r1].setdefault(r2, [])
                residues[r1][r2].append(viol)
                weights[r1][r2].append(w)
                                                       

                residues.setdefault(r2, {})
                weights.setdefault(r2, {})
                residues[r2].setdefault(r1, [])
                weights[r2].setdefault(r1, [])
                residues[r2][r1].append(viol)
                weights[r2][r1].append(w)
                
                    
        rms_per_resid = {}

        resid = residues.keys()
        resid.sort()

        matrix = []
        self.segids = []

        k = residues.keys()[0]
        m = k.getChain().getResidues()[-1].getNumber()
        self.rmat = zeros((m+1, m+1), float)
        
        for res1, res2s in residues.items():
            for res2, v in res2s.items():


                viol = array(v) 
                viol = power(viol * greater(viol, 0), 2) * array(weights[res1][res2])
                
                if sum(weights[res1][res2]) > 0:
                    rms  = sum(viol) / sum(weights[res1][res2])
                    rms = sqrt(rms)

                    # simple rules to distinct intra from inter
                    # if intra [a, b] with a > b
                    # if inter [a, b] with b > a
                    i1, i2 = self._get_order_from_type(res1,res2)
                    
                    s1, r1 = i1.getChain().getSegid(), i1.getNumber()
                    s2, r2 = i2.getChain().getSegid(), i2.getNumber()

                    matrix.append([r1,r2,rms])
                    self.segids.append([s1,s2])

                    self.rmat[r1][r2] = rms
                    self.rmat[r2][r1] = rms

        return resid, matrix

    # simple utiliy
    def _get_order_from_type(self, a, b):

        ra = a.getNumber()
        rb = b.getNumber()
        sa = a.getChain().getSegid()
        sb = b.getChain().getSegid()
        
        if sa <>  sb:
            
            if ra > rb:
                return b, a
            else:
                return a, b

        else:
            if ra > rb:
                return a, b
            else:
                return b, a
            
            
    def dumpRmsAnalysis(self, file = None):

        if file is None:
            file = self.text_report
            
        # Overall
        h = open(file, 'w')
        s = '''# Active restraints RMS violation
all restraints           %5.3f
ambiguous restraints     %5.3f
unambiguous restraints   %5.3f
''' % (self.rms, self.rms_a, self.rms_u)
        
        h.write(s)

        # profile
        s = '''
# RMS violation per residue
'''
        res = self.profile.keys()
        res.sort()
        for l in res:
            s += '%d    %5.3f\n' % (l, self.profile[l])
            
        h.write(s)


        # residue pairs
        s = '''
# RMS violation 2D map
'''
        
        for i in zip(self.rms_matrix, self.segids):
            s += '%3d %s - %3d %s    %5.3f\n' % (i[0][0], i[1][0], i[0][1], i[1][1], i[0][2])
                                                      
        h.write(s)
        h.close()


    def plotRmsProfile(self, file):

        import matplotlib.pylab as pylab
        
        
        pylab.axes([0.125, 0.7, 0.8, 0.2])
        x = [ r.getNumber() for r in self.mol.get_chains()[0].getResidues()]
        #x = self.profile.keys()
        #x.sort()
        #dat = [self.profile[i] for i in x]
        dat = zeros(len(x), float)
        for a in range(0, len(x)):
            dat[a] = self.profile.get(x[a]) or 0.
                  
        pylab.plot(x, dat, ls='-', c='b', marker="")
        pylab.grid(True)
        pylab.title('RMS violation residue profile')
        pylab.ylabel('RMS (A)')
        

        
    def plot2dRmsMap(self, file):

        import matplotlib.pylab as pylab
        from matplotlib import rcParams
	from numpy import ma
        rcParams['numerix'] = 'numpy'

        first = self.mol.get_chains()[0].getResidues()[0].getNumber()
        last = self.mol.get_chains()[0].getResidues()[-1].getNumber()
        
        self.rmat = self.rmat[first:,first:]
        X = ma.array(self.rmat, mask = equal(self.rmat, 0))
        pylab.axes([0.125, 0.1, 0.8, 0.5])
        #pylab.subplot(212)
        pylab.imshow(X, origin='lower', cmap=pylab.cm.Reds, interpolation='nearest', extent = (first,last, first, last))
        pylab.colorbar()
        pylab.title("RMS violation map")
        pylab.xlabel("Sequence")
        pylab.ylabel("Sequence")
        
        #pylab.savefig('/home/Bis/bardiaux/projects/relax/ph/toto.ps', papertype='a4')


    def plotPsReport(self, file = None):

        import matplotlib.pylab as pylab
        
        pylab.figure(figsize=(8,10))
        pylab.subplots_adjust(hspace=.9)
        
        if file is None:
            file = self.ps_report

        
            
        self.plotAllEps()

        it_n = self.iteration_n
        
        #pylab.subplots_adjust(hspace=0.5)
        pylab.figtext(0.3,0.95, 'RMS violation analysis for iteration %s' % str(it_n))
        pylab.savefig(file, papertype='a4')
        

    def plotAllEps(self):

        self.eps_profile = os.path.join(self.graphics_path, RMS_PS_PROFILE)
        self.eps_map = os.path.join(self.graphics_path,RMS_PS_MAP)
        self.eps_map_legend  = os.path.join(self.graphics_path, RMS_PS_MAP_LEGEND)
        
        self.plotRmsProfile(self.eps_profile)
        self.plot2dRmsMap(self.eps_map)

        

    def doRmsAnalysis(self):

        active = [p for p in self.peaks if p.isActive()]        
        inactive = [p for p in self.peaks if not p.isActive()]

        ambig = [p for p in active if p.isAmbiguous()]
        unambig = [p for p in active if not p.isAmbiguous()]

        self.mol = active[0].getContributions()[0].\
                   getSpinSystems()[0].getAtoms()[0].\
                   getResidue().getChain().getMolecule()

        # rms overall / ambig / unambig

        self.rms = self.getOverallRms(active)
        self.rms_u = self.getOverallRms(unambig)
        self.rms_a = self.getOverallRms(ambig)

        # rms profile

        self.profile = self.getRmsPerResidue(active)
        
        # rms 2D plot

        self.ids, self.rms_matrix = self.getRmsPerResiduePair(active)



    def go(self):

        try:
            # if something wrong there, no need to continue
            self.doRmsAnalysis()
        except:
            return

        # noe_restraints.rms
        self.dumpRmsAnalysis()

        # RmsPlot.ps
        if _use_matplotlib:
            self.plotPsReport()

if __name__ == '__main__':
    
    import sys, os
    path = os.environ['ARIA2']
    sys.path.insert(0, path + "/src/py/")
    
        
    from aria.tools import Load
    n = 8    
    out_path = '/home/Bis/bardiaux/devel/aria2.2_release/test/run2/structures/it' + str(n)
    
    file = out_path + '/noe_restraints.pickle'
    peaks = Load(file)

    rp = RmsReport(peaks, n, out_path, out_path)
    rp.go()
    

        

        
        

        

