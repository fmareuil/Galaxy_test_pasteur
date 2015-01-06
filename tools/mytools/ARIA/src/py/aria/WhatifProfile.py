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

CHECK_LIST = ['QUACHK', 'NQACHK', 'RAMCHK', 'C12CHK',
              'BBCCHK','HNDCHK', 'BMPCHK']

OTHER = ['BNDCHK', 'ANGCHK', 'PLNCHK', 'INOCHK']

LEGEND_LIST = ['1st generation packing quality Z-score (QUACHK)',
               '2nd generation packing quality Z-score (NQACHK)',
               'Ramachandran plot appearance Z-score (RAMCHK)',
               'Chi-1 chi-2 rotamer normality Z-score (C12CHK)',
               'Backbone conformation Z-score (BBCCHK)',
               'Bond lengths RMS Z-score (BNDCHK)',
               'Bond angles RMS Z-score (ANGCHK)',
               'Side chain planarity RMS Z-score (PLNCHK)',
               'Improper dihedral distribution RMS Z-score (HNDCHK)',
               'Inside/outside distribution RMS Z-score (INOCHK)',
               'Inter-atomic bumps (BMPCHK)']


FILENAME_WHATIF_PROFILE = 'whatif_profiles'
FILENAME_WHATIF_PROFILE_EPS = FILENAME_WHATIF_PROFILE + '.eps'
FILENAME_WHATIF_PROFILE_PS = FILENAME_WHATIF_PROFILE + '.ps'

_use_matplotlib = 0

try:
    import matplotlib
    _use_matplotlib = 1
    matplotlib.use('PS', warn=False)
except:
    pass


#FILENAME_WHATIF_PROFILE_TEX = FILENAME_WHATIF_PROFILE + '.tex'
    
class WhatifProfile:

    def __init__(self, wdir, whatIfExe = 'whatif'):
        
        self.checks = None
        self.data = None
        self.wdir = wdir
        self.is_whatcheck = 0

        if os.path.split(whatIfExe)[1] in ['whatcheck', 'DO_WHATCHECK.COM']:
            self.is_whatcheck = 1
            

    def parseWhatifChecks(self):

        from legacy.QualityChecks.QualityChecks import FILENAME_REPORT_WHATIF
        import re
        
        input = open(os.path.join(self.wdir, FILENAME_REPORT_WHATIF))
        
        fileID  = re.compile("^ID\s+\:\s+(\S+)")
        checkID = re.compile("^CheckID\s+\:\s+(\S+)")
        resid   = re.compile("\s+Name\s+\:\s(\S*)\-?\s*(\d+)\-(\w+)(\-\s*(CA))?")
        if self.is_whatcheck:
            resid   = re.compile("\s+Name\s+\:\s+(\d+)\s+;\s+\S*\s+;\s+(\d+)\s+;\s+(\w+)\s+;\s+\S*\s+(;\s+CA)?")        
        value   = re.compile("\s+Value\s+\:\s+(\-?\d+\.\d+)")

        checks = {}

        self.check_list = []
        # parsing
        while 1:

            line = input.readline()

            fID = fileID.search(line)
            cID = checkID.search(line)
            res = resid.search(line)    
            val = value.search(line)

            if fID:

                file = fID.group(1)
                checks.setdefault(file, {})

            if cID:

                ch = cID.group(1)
                checks[file].setdefault(ch, {})
                if ch not in self.check_list and ch in CHECK_LIST:
                    self.check_list.append(ch)

            if res:
                resn = res.group(2)

            if val:
                v = float(val.group(1))
                try:
                    x = checks[file][ch][resn]
                    continue
                except:
                    checks[file][ch][resn] = v

            if not line:
                break

        self.checks = checks
        


    def getProfilesData(self):
        
        import numpy

        self.plength = len(self.checks[self.checks.keys()[0]]['NQACHK'])
        self.pids = map(lambda a: int(a), self.checks[self.checks.keys()[0]]['NQACHK'].keys())
        self.pids.sort()
        
        #print self.plength

        scores = {}
        extrem = {}
        sd = {}

        for elt in self.check_list:

            score = numpy.array([ cv for k,v in self.checks.items() for ck, cv in v.items() if ck == elt])

            new_score = []

            #print elt
            for s in score:
                struc_score = numpy.zeros(self.plength, numpy.float)
                
                for g in range(self.plength):
                    try:
                        #struc_score[g] = s[str(g+1)]
                        struc_score[g] = s[str(self.pids[g])]
                        
                    except:
                        pass
                new_score.append(struc_score)
                
            new_score = numpy.array(new_score)
            #if new_score
            #print new_score
            scores[elt] = numpy.mean(new_score, 0)
            extrem[elt] = [numpy.min(new_score, 0), numpy.max(new_score, 0)]
            if len(new_score) > 1:
                sd[elt] =  numpy.std(new_score, 0)
            else:
                sd[elt] = numpy.zeros(len(new_score[0]), numpy.float)
                           

        self.data = [scores, extrem, sd] 


    def plotProfiles(self, output_ps):
        
        #import biggles
        import matplotlib.pylab as pylab
        
        from matplotlib import rcParams
        #from numpy import ma
        rcParams['numerix'] = 'numpy'            

        from math import ceil

        # Test a one-col table
        n_l = len(self.check_list)

        f = pylab.figure(figsize=(8,10))
        pylab.figtext(0.3,0.95,"Whatif scores residue profiles")
        for e in range(len(self.check_list)):
            
            name = self.check_list[e]
            if e >= n_l:
                x = e - n_l
                y = 1
            else:
                x = e
                y = 0

            # test suite
            y = 0
            x = e
            scores = self.data[0]
            pylab.subplot(len(self.check_list), 1, e+1)
            pp = self.unitaryPlot(name, scores[name])
 
        pylab.subplots_adjust(hspace=0.5)
        #t1.write_eps(output_ps)
        pylab.savefig(output_ps, papertype='a4')
        

    def unitaryPlot(self, t, data):

        #import biggles
        import matplotlib.pylab as pylab

        linet = ['dotted', 'solid']
        profile = []


        if type(data) == type([]):
            profile = data        
        else:
            profile.append(data)

        #x = range(1, len(profile[0])+1)
        x = self.pids

        n = 1
        for d in profile:


            pylab.plot(x, d, ls='-', c='b', marker=".")
            pylab.grid(True)
            pylab.ylabel(t)
            pylab.title(t)
            pylab.xlim(x[0], x[-1])
            
            n+=1



    def writeTextProfiles(self, dst):
        """
        write profiles in text file
        resid         check1             check2
                  avg max min sd     avg max min sd
        """
        
        line = """# Whatif scores residues profiles

# Check  : name of the whatif analysis
# Res    : residue number
# Mean   : mean score values among all conformers
# Sd     : unbiased standard deviation
# Max    : maximum score
# Min    : minimum score

"""
        
        
        line += '#%6s\t%3s\t%6s\t%6s\t%6s\t%6s\n\n' % ('Check', 'Res', 'Mean', 'Sd', 'Min', 'Max')
        
        for e in self.check_list:
            for r in range(self.plength):

                idn = self.pids[r]
                
                line += '%6s\t%3d\t%6.3f\t%6.3f\t%6.3f\t%6.3f\n' % (e,idn,self.data[0][e][r],self.data[2][e][r], self.data[1][e][0][r], self.data[1][e][1][r])
            line += '\n'

        f = open(dst, 'w')
        f.write(line)
        f.close()
        

    def makeProfiles(self):

        # Parse whatif output residue file
        self.parseWhatifChecks()

        # average scores on all        
        self.getProfilesData()

        #try:
        # plot all checks profiles
        if _use_matplotlib:
            output_ps = os.path.join(self.wdir, FILENAME_WHATIF_PROFILE_PS)
            self.plotProfiles(output_ps)

            # if we are not on water refinement, change location of ps and eps to /graphics
            from aria.tools import copy_file
            src = output_ps
            dst = os.path.join(os.path.join(self.wdir, 'graphics'), FILENAME_WHATIF_PROFILE_PS)
            if os.path.exists(os.path.join(self.wdir, 'graphics')):
                try:
                    copy_file(src, dst)
                    os.remove(src)
                except IOError:
                    s = 'File "%s" does not exist or cannot be accessed.'
                    print IOError, s % src
            
        #except:
        #    print "Plotting abording. Check if module biggles/ps installed."

        # write ASCII whatif residue profiles 
        output = os.path.join(self.wdir, FILENAME_WHATIF_PROFILE)
        self.writeTextProfiles(output)


    


if __name__ == '__main__':

    import sys, os

    path = os.environ['ARIA2']

    sys.path.insert(0, path+'/src/py')

    p = WhatifProfile(sys.argv[1])
    p.makeProfiles()
