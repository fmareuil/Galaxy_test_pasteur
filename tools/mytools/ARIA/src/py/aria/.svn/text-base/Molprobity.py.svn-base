"""
Authors: Bardiaux Benjamin
         Structural Biology Unit, FMP, Berlin
        
         Copyright (C) Benjamin Bardiaux
         No warranty implied or expressed.
         All rights reserved.

$Author: bardiaux $
$Revision: 1.1.1.1 $
$Date: 2010/03/23 15:27:24 $
"""

import os

class MolprobityClashlist:


    def __init__(self, clashlistExe):

        self.clashlistExe = os.path.expanduser(clashlistExe)
        
    def runClashlist(self, pdbfiles, outputFile):

        cmd_line = "%s %s"

        outF = open(outputFile, 'w')
        
        for pdbfile in pdbfiles:
            outF.write("#" + pdbfile + "\n")
            cmd = cmd_line % (self.clashlistExe, pdbfile)
            output = os.popen(cmd).read()
            outF.write(output)

        outF.close()

        scores = self.readClashlist(outputFile)

        return scores
             

    def readClashlist(self, inputFile):

        scores = []
        inF = open(inputFile)
        for line in inF:
            if not line:
                break
            elif line.startswith("#sum2"):
                score = float(line.split(":")[2].split('clashscore')[0])
                scores.append(score)
        inF.close()

        return scores
                           

    def doStatistics(self, scores):

        db = os.path.dirname(self.clashlistExe)
        db = os.path.dirname(db)
        db = os.path.join(db, "lib/clashscore.db.tab")

        dbScores = self.readClashlistDb(db)

        # log-norm model for Z-scores
        from numpy import log, mean, std, array, compress, greater, greater_equal, fix

        dbScoresLog = log(compress(greater(dbScores, 0.), dbScores))
        mu = mean(dbScoresLog)
        sd = std(dbScoresLog)

        Z = (mu - log(scores))/sd

        # percentile rank
        Ndb = len(dbScores)

        pc = array([sum(greater_equal(dbScores, sc)) for sc in scores])
        pc = pc*100./Ndb
        pc = fix(pc)

        return Z.tolist(), pc.tolist()

    def readClashlistDb(self, db):

        dbScores= []
        dbF = open(db)
        for line in dbF:
            if not line:
                break
            elif not line.startswith("#"):
                d = line.split(":")
                score = float(d[2])
                dbScores.append(score)

        dbF.close()
        
        return dbScores
            


if __name__ == '__main__':


    import sys
    
    exe = '~/Progs/molprobity3/bin/clashlist'

    mpcl = MolprobityClashlist(exe)

    pdbs = sys.argv[1:]
    outF = 'molpro.dat'
    
    scores = mpcl.runClashlist(pdbs, outF)

    Z, pc = mpcl.doStatistics(scores)

    print scores
    print Z
    print pc
