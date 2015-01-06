"""

          ARIA -- Ambiguous Restraints for Iterative Assignment

                 A software for automated NOE assignment

                               Version 2.3


Copyright (C) Benjamin Bardiaux, Michael Habeck, Therese Malliavin,
              Wolfgang Rieping, and Michael Nilges

All rights reserved.


NO WARRANTY. This software package is provided 'as is' without warranty of
any kind, expressed or implied, including, but not limited to the implied
warranties of merchantability and fitness for a particular purpose or
a warranty of non-infringement.

Distribution of substantively modified versions of this module is
prohibited without the explicit permission of the copyright holders.

$Author: bardiaux $
$Revision: 1.1.1.1 $
$Date: 2010/03/23 15:27:24 $
"""





## Settings
from aria.Settings import Settings

########## BARDIAUX
class CmapSettings(Settings):

    def create(self):

        from aria.TypeChecking import DICT, LIST

        from aria.Settings import TypeEntity, ChoiceEntity, NonNegativeInt, Path, AbsolutePath, String

        kk = {}
        entity = TypeEntity(DICT)
        entity.set({})
        kk['xy'] = entity
        entity = TypeEntity(LIST)
        entity.set([])
        kk['posit'] = entity
        entity = TypeEntity(DICT)
        entity.set({})
        kk['plist'] = entity
        entity = TypeEntity(DICT)
        entity.set({})
        kk['clist'] = entity
        
        kk['file'] = AbsolutePath()
        kk['it'] = NonNegativeInt()

        kk['selection'] = String()
        
        return kk

    def create_default_values(self):

        d = {}
        d['xy'] = {}
        d['posit'] = []
        d['plist'] = {}
        d['clist'] = {}
        d['file'] = ""
        d['selection'] = ""
        
        return d
    
class ContribSettings(Settings):


    def create(self):
        
        from aria.Settings import NonNegativeInt, String
        d = {}
        
        d['ariapeak'] = NonNegativeInt()
        d['id'] = NonNegativeInt()
        d['dist'] = String()
        d['weight'] =  String()
        d['res1'] = NonNegativeInt()
        d['at1'] = String()        
        d['seg1'] = String()
        d['res2'] = NonNegativeInt()
        d['at2'] = String()        
        d['seg2'] = String()        
        
        return d

    def create_default_values(self):

        d = {}

        d['ariapeak'] = 0        
        d['id'] = 0
        d['dist'] = '%5.3f' % 0.
        d['weight'] = '%5.3f' % 0.
        d['res1'] = 0
        d['at1'] = ""       
        d['seg1'] = ""
        d['res2'] = 0
        d['at2'] = ""     
        d['seg2'] = ""      
        d['seg2'] = "" 

        return d

class PeakSettings(Settings):

    def create(self):
        
        from aria.Settings import NonNegativeInt, String, ChoiceEntity
        from aria.ariabase import YES, NO
        d = {}
        
        d['id'] = NonNegativeInt()
        d['ref_peak'] = NonNegativeInt()
        d['spec'] = String()

        val = ('active','inactive')
        d['state'] = ChoiceEntity(val)        
        
        d['dist'] = String()
        d['up'] = String()
        d['low'] = String()
        
        d['weight'] = String()
                
        d['deg_viol'] = String()
        d['dist_avg'] = String()

        val = (YES, NO)
        d['viol'] = ChoiceEntity(val)
        
        val = ('unambiguous','ambiguous')        
        d['type'] = ChoiceEntity(val)
        
        return d

    def create_default_values(self):

        d = {}
        d['id'] = 0
        d['ref_peak'] = 0
        d['spec'] = ""

        d['state'] = 'inactive'        
        
        d['dist'] = '%5.3f' % 0.
        d['up'] = '%5.3f' % 0.
        d['low'] = '%5.3f' % 0.
        
        d['weight'] = '%5.3f' % 0.
                
        d['deg_viol'] = '%5.2f %%' % 0.
        d['dist_avg'] = '%5.3f' % 0.
        
        d['viol'] = 'no' 
    

        return d

class PeakListSettings(Settings):




    def create(self):
        
        from aria.TypeChecking import LIST
        from aria.Settings import NonNegativeInt, TypeEntity
        d = {}

        
        d['id'] = NonNegativeInt()

        entity = TypeEntity(LIST)
        entity.set([])
        d['pk_list'] = entity
        
        return d

    def create_default_values(self):

        d = {}
        
        d['id'] = 0
        d['pk_list'] = []   

        return d

class ContribListSettings(Settings):



    def create(self):
        
        from aria.TypeChecking import LIST
        from aria.Settings import NonNegativeInt, TypeEntity
        d = {}

        
        d['id'] = NonNegativeInt()

        entity = TypeEntity(LIST)
        entity.set([])
        d['c_list'] = entity
        
        return d

    def create_default_values(self):

        d = {}
        
        d['id'] = 0
        d['c_list'] = []   

        return d
    
class DisplaySettings(Settings):



    def create(self):
        from aria.Settings import String
        d = {}

        
        d['display'] = String()

        
        return d

    def create_default_values(self):

        d = {}
        
        d['display'] = 'all'   

        return d

######## PEAK LIST CALUCLATION AND CMAP DETREMINATION
def get_cmap_list(ap_list, sel):
    
    if sel == 'all':
        peaks = [p for p in ap_list if p.isActive()]
    elif sel == 'ambig':
        peaks = [p for p in ap_list if p.isActive() and p.isAmbiguous()]
    elif sel == 'unambig':
        peaks = [p for p in ap_list if p.isActive() and not p.isAmbiguous()]        
        
    
    n=0
    cmap = []
    pk_list = {}
    con_list = {}
   
    for p in peaks:

        contr = [c for c in p.getContributions() if c.getWeight() > 0.]

        for c in contr:

            cple = []
            
            for spin_sys in c.getSpinSystems():
                at = spin_sys.getAtoms()[0]
                re = at.getResidue().getNumber()
                cple.append(re)

            cples = [cple]
            if cple[0] <> cple[1]:
                cples.append([cple[1],cple[0]])

            for cple in cples:

                if cple not in cmap:
                    cmap.append(cple)
                    con_list.setdefault(n,[setContrib([c,p.getId()])])
                    pk_list.setdefault(n,[setPeak(p)])
                    n+=1
                else:
                    i = cmap.index(cple)
                    if [c,p.getId()] not in con_list[i]:
                        con_list[i].append(setContrib([c,p.getId()]))
                    if setPeak(p) not in pk_list[i]:
                        pk_list[i].append(setPeak(p))                
                    
                
##     uniq=reduce(lambda l, x: x not in l and l.append(x) or l, cmap , [])

    return cmap, pk_list, con_list
        
def setPeak(p):

##     from gui_beta import PeakSettings

    from aria.ariabase import YES, NO

    YN = [NO, YES]
    AC = ['inactive','active']
    AMB = ['unambiguous','ambiguous']
    
    pp = {}
    
    pp['id'] = p.getId()
    pp['ref_peak'] = p.getReferencePeak().getNumber()
    pp['spec'] = p.getReferencePeak().getSpectrum().getName()
    pp['state'] = AC[p.isActive()]
    pp['dist'] = '%5.3f' % p.getDistance()
    pp['up'] = '%5.3f' % p.getUpperBound()
    pp['low'] = '%5.3f' % p.getLowerBound()
    pp['weight'] = '%5.3f' % p.getWeight()
    pp['deg_viol'] = '%5.2f %%' % (p.analysis.getDegreeOfViolation() * 100)
    pp['dist_avg'] = '%5.3f' % p.analysis.getAverageDistance().getValue()
    pp['viol'] = YN[p.analysis.isViolated()]
    pp['type'] = AMB[p.isAmbiguous()]

    return pp

def setContrib(c):

##     from gui_beta import PeakSettings
        
    
    pp = {}
    
    cc = c[0]
    
    INT = ['INTER','INTRA']
    
    pp['id'] = cc.getId()
    pp['ariapeak'] = c[1]
    if cc.getAverageDistance().getValue() is not None:
        pp['dist'] = '%5.3f' % cc.getAverageDistance().getValue()
    else:
        pp['dist'] = 'None'
    pp['weight'] = '%5.3f' % cc.getWeight()
    
    at1 = cc.getSpinSystems()[0].getAtoms()
    atoms1 = pseudoString([a.getName() for a in at1])
    
    
    at2 = cc.getSpinSystems()[1].getAtoms()
    atoms2 = pseudoString([a.getName() for a in at2])
    
    pp['res1'] =at1[0].getResidue().getNumber()
    pp['at1'] = ' / '.join(atoms1)
    pp['seg1'] = at1[0].getSegid()
    pp['res2'] = at2[0].getResidue().getNumber()
    pp['at2'] = ' / '.join(atoms2)
    pp['seg2'] = at2[0].getSegid()
    
    return pp

def pseudoString(at):

    import re
    seen = []
    s = '(\D+)\d+'
    m = re.compile(s)
    for a in at:
        x = m.search(a)
        if x:
            _a = x.group(1) + '*'
            if _a not in seen:
                seen.append(_a)
        else:
            if a not in seen:
                seen.append(a)            

    return seen


class Exporter:

    def __init__(self):

        pass

    def write(self, filename):

        pass

    

    
