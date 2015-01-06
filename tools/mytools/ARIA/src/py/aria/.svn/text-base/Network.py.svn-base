"""
Authors: Bardiaux Benjamin
         Institut Pasteur, Paris
         IBPC, Paris
        
         Copyright (C) 2005 Michael Habeck,
         Wolfgang Rieping and Benjamin Bardiaux
         No warranty implied or expressed.
         All rights reserved.

$Author: bardiaux $
$Revision: 1.1.1.1 $
$Date: 2010/03/23 15:27:24 $
"""

from aria.ariabase import *
from aria.Settings import Settings
from aria.xmlutils import XMLElement, XMLBasePickler
import aria.TypeChecking as TCheck
from aria.Chain import TYPE_NONPOLYMER

import numpy

from time import clock

from aria.AriaPeak import TextPickler
from aria.AriaPeak import ASSIGNMENT_TYPE_DICT, NA, \
     HEADER_PROJECT, HEADER_ASSIGNMENT_TYPE, \
     HEADER_SEQUENCE_SEPARATION, HEADER_RESTRAINT_DEFINITION, \
     HEADER_RESTRAINT_ACTIVE



HEADER_SEQUENCE_SEPARATION = \
"""
# sep:      sequence separation s: I:      s == 0  (intra-residual)
#                                  Q:      s == 1  (sequential)        
#                                  S: 2 <= s <= 3  (short)
#                                  M: 4 <= s <= 5  (medium)
#                                  L:      s >  5  (long)
#                                  i: inter-monomer
"""[1:-1]


HEADER_DICT = {'project': HEADER_PROJECT,
               'assignment_type': HEADER_ASSIGNMENT_TYPE,
               'sequence_separation': HEADER_SEQUENCE_SEPARATION,
               'restraint_definition': HEADER_RESTRAINT_DEFINITION,
               'restraint_active': HEADER_RESTRAINT_ACTIVE}

HEADER_ABBREVIATIONS = \
("""
#
# Abbreviations:
#
%(restraint_definition)s
%(restraint_active)s
#
#
%(assignment_type)s
#
""" % HEADER_DICT)[1:-1]


HEADER_ALL = \
"""
#
# List of distance restraints.
#
# Created by Aria 2.3, %(creation_date)s
#
%(project)s
#
# Restraints used during calculation: %(n_active)d
# Violated: %(n_violated)d
#
%(abbreviations)s
%(sequence_separation)s
#
# n_c:      The number of contributions. (see noe_restraints.assignments for
#           explicit list of contributions).
#
# net_res:  Network-anchoring score per residue.
#
# net_ato:  Network-anchoring score per atom.
#
"""[1:]


class NetworkScoreTextPickler(TextPickler):

    def encode_common(self, ap):
        
        distance_format = '%.2f'

        number = '%d' % ap.getId()
        
        rp = ap.getReferencePeak()

        x = rp.getNumber()

        try:
            ref_peak_number = '%d' % x
        except:
            ref_peak_number = NA

        x = rp.getSpectrum().getName()

        try:
            ref_peak_spectrum = str(x)
        except:
            ref_peak_spectrum = NA

        x = ap.isActive()

        if x:
            active = YES
        else:
            active = NO

        at = rp.getAssignmentType()
        assignment_type = ASSIGNMENT_TYPE_DICT[at]

        # BARDIAUX
        net = ap._network
        net_res = '%.2f' % ap._network['residue']
        net_ato = '%.2f' % ap._network['atom']

        values = ref_peak_spectrum, ref_peak_number, number, \
                 active,  net_res, net_ato, assignment_type

        return list(values)

    def encode(self, ap):

        values = self.encode_common(ap)

        ## contributions

        contributions = ap.getContributions()

        ## take only active contributions

        contributions = ap.getActiveContributions()

        if len(contributions) == 1:

            ## get sequence separation
            ## in case of multuple spin-pairs,
            ## we just take the first one, since all are
            ## involve the same two residues

            atom1, atom2 = contributions[0].getSpinPairs()[0].getAtoms()

            if atom1.getSegid() <> atom2.getSegid():
                # we have an inter
                values.append('1') # n_c
                values.append('i')
                return values
                 
            seq_pos1 = atom1.getResidue().getNumber()
            seq_pos2 = atom2.getResidue().getNumber()

            seq_sep = abs(seq_pos1 - seq_pos2)

            ## intra-residue

            if seq_sep == 0:
                descr = 'I'

            ## sequential
                
            elif seq_sep == 1:
                descr = 'Q'

            ## TODO: are these the correct values?

            ## short range

            elif seq_sep <= 3:
                descr = 'S'

            ## medium range

            elif seq_sep <= 5:
                descr = 'M'

            else:
                descr = 'L'

            values.append('1') # n_c
            values.append(descr)
            

        ## multiple contributions

        else:
            values.append(str(len(contributions)))
            values.append('-') # sep

        return values

    def dumps(self, ap):
        return '\n'.join(self.encode(ap))
    
    
class NetworkAnchoringTextPickler(TextPickler):

    HEADER_COMMON = ['ref_spec', 'ref_no', 'id', 'active', 'net_res', 'net_ato', 'a_type']

    COLUMNS = {'all' : HEADER_COMMON + ['n_c', 'sep'],}

    HEADER = {'all' : HEADER_ALL,}

    def __init__(self, settings):
        #check_type(settings, 'AriaPeakListTextPicklerSettings')
        TextPickler.__init__(self, settings = settings)
    
    def get_column_header(self, _type):
        """
        _type is 'ambig' or 'unambig'
        """

        if not _type in ('ambig', 'unambig', 'all'):
            s = 'Header for peak-type "%s" not known.' % _type
            self.error(TypeError, s)

        return list(self.COLUMNS[_type])

    def encode(self, peak_list, header):

        pickler = NetworkScoreTextPickler()
        all = map(pickler.encode, peak_list)
    
        ## add header

        if not len(all):
            return header

        if len(header) <> len(all[0]):
            s = 'Number of columns must match header-length.'
            self.error(Exception, s)

        header[0] = '# ' + header[0]
        
        ## show additional information

        active = [p for p in peak_list if p.isActive()]

        n_violated = len([p for p in active if p.analysis.isViolated()])

        d = self._compile_header_dict()
        
        d['n_violated'] = n_violated
        d['n_active'] = len(active)
        d['abbreviations'] = HEADER_ABBREVIATIONS
        text = self.format_output(all, header = header)

        ## add \n
        text = [line + '\n' for line in text]

        ## make string

        text = ''.join(text)

        return text, d
    
    def _write(self, s, filename, gzip = 0):
        import os
        
        if s is None:
            import aria.tools as tools
            
            tools.touch(filename)
            return

        if gzip:
            from aria.tools import gzip_open as open_func
        else:
            open_func = open

        filename = os.path.expanduser(filename)

        f = open_func(filename, 'w')
        f.write(s)
        f.close()

    def _compile_header_dict(self):
        from aria.Singleton import ProjectSingleton
        import time
        from copy import copy
        
        project = ProjectSingleton()
        project_settings = project.getSettings()

        infra = project.getInfrastructure()
        run_path = infra.get_run_path()

        d = {'date': project_settings['date'],
             'project': project_settings['name'],
             'run': project_settings['run'],
             'author': project_settings['author'],
             'working_directory': run_path}

        x = copy(HEADER_DICT)
        x['project'] %= d
        x['creation_date'] =time.ctime()

        return x

    def dump_network(self, peak_list, filename, gzip = 0):
        
        if peak_list:

            header = self.get_column_header('all')
            text, d = self.encode(peak_list, header)
            d.update(self._compile_header_dict())
                     
            header = (self.HEADER['all'] % d)[1:]
            s = header + text
#            s = header.replace('\n\n','\n') + text

        else:
            s = None
        
        return self._write(s, filename, gzip)

class NetworkPsPickler:
        
    def __init__(self, network):

        self.peaks = network.peaks
        self.p_id = network._protons_id
        self.net_res = network.residue_score
        self.mol = network.molecule
        self.it_n = network.iteration.getNumber()


    def get_matrix(self):

        # since we just support symmetric dimer
        n_chains = len(self.mol.get_chains())

        #max_res = len([r for c in self.mol.get_chains() for r in c.getResidues()])
        max_res = [c.getResidues()[-1].getNumber() for c in self.mol.get_chains() \
                   if c.getType() != TYPE_NONPOLYMER]

        from aria.Singleton import ProjectSingleton
        from aria.DataContainer import DATA_SYMMETRY
        project = ProjectSingleton()
        sym_settings = project.getData(DATA_SYMMETRY)[0]        

        if n_chains < 2 or (n_chains > 1 and sym_settings['symmetry_type'] not in ["C2","C3","D2","C5"]):
            # monomeric prot or hetero dimer
            
            matrix = numpy.zeros((max_res[0]+1, max_res[0]+1), numpy.float)

            for k, r_net in self.net_res.items():
                r1, r2 = map(lambda a: a.getNumber(), k)
                matrix[r1,r2] = r_net
                matrix[r2,r1] = r_net

            return matrix, None

        else:
            # homo-dimer
            matrix_a = numpy.zeros((max_res[0]+1, max_res[0]+1), numpy.float)
            matrix_r = numpy.zeros((max_res[0]+1, max_res[1]+1), numpy.float)
            
            for k, r_net in self.net_res.items():
                r1, r2 = map(lambda a: a.getNumber(), k)
                s1, s2 = map(lambda a: a.getChain().getSegid(), k)

                if s1 <> s2:
                    matrix_r[r1,r2] = r_net
                    matrix_r[r2,r1] = r_net
                else:
                    matrix_a[r1,r2] = r_net
                    matrix_a[r2,r1] = r_net

            return matrix_a, matrix_r
                    

    def plot_matrix(self):
        
        # mask zero-values
        from matplotlib import rcParams
        from numpy import ma
        rcParams['numerix'] = 'numpy'
            
        pylab = self.pylab
        msg = ""
        
        matrix_a, matrix_r = self.get_matrix()

        first_res = [c.getResidues()[0].getNumber() for c in self.mol.get_chains() if c.getType() != TYPE_NONPOLYMER]
        max_res = [c.getResidues()[-1].getNumber() for c in self.mol.get_chains() if c.getType() != TYPE_NONPOLYMER]
        
        if matrix_r is not None:
            
            ax1 = pylab.subplot(2,1,1)

            #matrix = matrix_r[1:,1:]
            matrix = matrix_r[first_res[0]:,first_res[1]:]
            
            X = ma.array(matrix, mask = numpy.equal(matrix, 0.))

            xyticks = (first_res[0], max_res[0], first_res[1], max_res[1])
            
            kw = {'origin':'lower',
                  'interpolation':'nearest',
                  'aspect' :  'equal',
                  'extent' : xyticks}
            
            pylab.imshow(X, cmap=pylab.cm.Reds, **kw)

            pylab.grid()
            pylab.colorbar(orientation = 'vertical')
            pylab.ylabel("Residue Number (Inter-molecular)")
            #pylab.setp( ax1.get_xticklabels(), visible=False)
            
            pylab.subplot(212)#, sharex=ax1)

            #pos = pylab.axes([0.85, 0.1, 0.04, 0.8])
            #pylab.colorbar(cax = pos)#, orientation = 'horizontal')

            msg = " (Intra-molecular)"

        matrix = matrix_a[first_res[0]:,first_res[0]:]
        #matrix = matrix_a[1:,1:]
        
        X = ma.array(matrix, mask = numpy.equal(matrix, 0.))
        
        xyticks = (first_res[0], max_res[0], first_res[0], max_res[0])
            
        kw = {'origin':'lower',
              'interpolation':'nearest',
              'aspect' :  'equal',
              'extent' : xyticks}
        
        pylab.imshow(X, cmap=pylab.cm.Reds, **kw)

        
        if len(msg):
            orientation = 'vertical'
        else:
            orientation = 'horizontal'
        pylab.colorbar(orientation = orientation)
        pylab.grid()
        pylab.xlabel("Residue Number")
        pylab.ylabel("Residue Number" + msg)

            

    def plot_profile(self, type, n):

        pylab = self.pylab

        if type not in ['residue', 'atom']:
            return

        colors = {'residue' : 'b',
                  'atom' : 'r'}

        scores = [p._network[type] for p in self.peaks]
        nbins  = int(max(scores))

        #nbins = 1 + int(numpy.log(len(scores))/numpy.log(2))
        nbins = int(1.0 + 3.3 * numpy.log(len(scores)))

        pylab.subplot(2, 1, n)
        pylab.hist(scores, bins = nbins +1, facecolor = colors[type])
        pylab.xlabel("Network Anchoring score per %s" % type)
        pylab.ylabel("Number of Peaks")

    

    def plot(self, path):

        try:
            import matplotlib
            matplotlib.use('PS', warn=False)
        except:
            return

        import matplotlib.pylab as pylab
        
        self.pylab = pylab
        
        pylab.figure(num=1, figsize=(8,11))
        pylab.clf()
        pylab.figtext(0.3,0.95, 'Network Anchoring for iteration %s' % str(self.it_n))
        pylab.figtext(0.3,0.90, 'Network Anchoring scores distribution')
        
        self.plot_profile('residue', 1)
        self.plot_profile('atom', 2)
        pylab.subplots_adjust(top = 0.85)
        
        pylab.figure(num=2, figsize=(8,11))
        pylab.clf()
        pylab.figtext(0.3,0.95, 'Residue-wise Network Anchoring scores for iteration %d' % self.it_n)
        self.plot_matrix()

        pylab.figure(1)
        pylab.savefig(path +'_dist.ps', papertype='a4', dpi = 72)
        pylab.figure(2)
        pylab.savefig(path + '_2D.ps', papertype='a4', dpi = 72)
        
        
class NetworkSettings(Settings):

    def create(self):

        from aria.Settings import NonNegativeFloat
        from aria.Settings import YesNoChoice

        d = {}

        # public settings
        descr = "Network anchoring removes restraints which are not surrounded by a network of active restraints."        
        d['enabled'] = YesNoChoice(description = descr)
        
        descr = "High network-anchoring score per residue for a peak to be active."        
        d['high_residue_threshold'] = NonNegativeFloat(description = descr)
        
        descr = """Minimal network-anchoring score per residue for a peak to be active. (In combination with \"min_atom_threshold\")"""        
        d['min_residue_threshold'] = NonNegativeFloat(description = descr)        

        descr = """Minimal network-anchoring score per atoms for a peak to be active. (In combination with \"min_residue_threshold\")"""        
        d['min_atom_threshold'] = NonNegativeFloat(description = descr)
        
        # private
        descr = "Maximal distance for covalent inter-proton distance."        
        d['distance_max'] = NonNegativeFloat(description = descr)

        descr = "Maximal network anchoring score for covalent distance."        
        d['v_max'] = NonNegativeFloat(description = descr)
        
        descr = "Minimal network anchoring score for intraresidual/sequential distance."        
        d['v_min'] = NonNegativeFloat(description = descr)
        
        return d

    def create_default_values(self):
        d = {}

        d['enabled'] = NO
        d['high_residue_threshold'] = 4.
        d['min_residue_threshold'] = 1.0
        d['min_atom_threshold'] = 0.25
        
        d['distance_max'] = 5.5
        d['v_max'] = 1.0
        d['v_min'] = 0.1
        
        return d
    
class CovalentConstraint:

    def __init__(self, id, atom1, atom2, distance):

        self.atom1 = atom1
        self.atom2 = atom2
        self.distance = distance
        self.id = id

    def getId(self):
        return self.id
    
    def getScore(self):
        return 0.

    def getAtoms(self):
        return (self.atom1, self.atom2)

    def getDistance(self):
        return self.distance

    def __str__(self):

        s = "CovalentConstraint(id=%d, atoms=%s, d=%5.3f)" % (self.id, self.getAtoms(), self.distance)
        return s
        
class NetworkAnchoring(AriaBaseClass):
        
    def __init__(self, settings):
        
        TCheck.check_type(settings, 'NetworkSettings')

        AriaBaseClass.__init__(self)

        self.setSettings(settings)

        self.anchoring = None
        self.peaks = None

        self.getSettings()['v_min'] = 0.1
        self.getSettings()['v_max'] = 1.0        
        self.getSettings()['distance_max'] = 5.5
        
        
    def setup(self):
        """
        Setup some lists and matrices.
        """

        from sets import Set

        if self.anchoring is not None:
            
            # if we already have a network, just recreate self._c_id with copied contribuitions
            self.message('Retrieving Network ...')
            
            self._c_id = {}
            self._c_id[-1] = [] # covalent
            
            for p in self.peaks:
                for c in p.getContributions():
                    for sp in c.getSpinPairs():
                        sid = sp.getId() + 1
                        self._c_id.setdefault(sid, Set())
                        self._c_id[sid].add(c)

            
            self.addDistanceRestraints()
            
            return 1
                
            
        # if we run network_anchoring for 1st time, create all list and spinpair matrices
        self.message('Initializing ...')

        if not self.peaks:
            return 0
        
        # list with all protons
        if self._is_noesy_only:
            self._protons_id = [a for c in self.molecule.get_chains() for r in c.getResidues() \
                                for a in r.getAtoms() if a.isProton()]
        else:
            self._protons_id = [a for c in self.molecule.get_chains() for r in c.getResidues() \
                                for a in r.getAtoms() if a.isProton() or a.getType() in ['N','C']]
        
        self._protons_id.sort(lambda a,b: cmp(a.getId(), b.getId()))

        # dict with protons id as key, and indices in self._protons_id as values
        self._protons_num = {}
        for a in range(0, len(self._protons_id)):
            self._protons_num[self._protons_id[a].getId()] = a

        # list with protons residues number
        # add chain levels to residues numbering
        self._residues_num = {}
        for c in self.molecule.get_chains():
            cid = c.getSegid()
            self._residues_num[cid] = [a.getResidue().getNumber() for a in self._protons_id]# if a.getSegid() == cid]

        # dict with residues number as key and list of protons ids as values
        self._residues_id = {}

        for c in self.molecule.get_chains():
            cid = c.getSegid()
            self._residues_id[cid] = {}
            
        for a in range(0, len(self._protons_id)):
            r, cid = self._protons_id[a].getResidue().getNumber(), self._protons_id[a].getSegid()
            self._residues_id[cid].setdefault(r, [])
            self._residues_id[cid][r].append(a)


        # dict with SpinPair.getId() + 1 as key and Set of contributions as values
        self._c_id = {}
        self._c_id[-1] = []
        
        # dict with SpinPair.getId() + 1 as key and spinpair as values
        self.spinpairs = {}
        
        for p in self.peaks:
            for c in p.getContributions():
                for sp in c.getSpinPairs():
                    sid = sp.getId() + 1
                    self._c_id.setdefault(sid, Set())
                    self._c_id[sid].add(c)

                    if not self.spinpairs.has_key(sid):
                        self.spinpairs[sid] = sp


        # add additional distance restraints
        self.addDistanceRestraints()
        
        # matrix to hold wether 2 protons are connected with spinpair(1), covalent(2) or not connected(0)
        self._sp = numpy.zeros((len(self._protons_id), len(self._protons_id)))
        
        # matrix to store the id of the spinpair connecting 2 atoms
        self._sp_id = numpy.zeros((len(self._protons_id), len(self._protons_id)))

        # matrix to store covalent score of a spinpair
        self._sp_cov_scores = numpy.zeros((len(self._protons_id), len(self._protons_id)))

        # matrix to store sum of contributions volumes of each spinpair
        self._sp_sum_scores = numpy.zeros(len(self.spinpairs.keys()) , numpy.float)
        
        for spid, sp in self.spinpairs.items():
            a, b = sp.getAtoms()
            a, b = self._protons_num[a.getId()], self._protons_num[b.getId()]

            self._sp[a][b] = 1
            self._sp[b][a] = 1

            self._sp_id[a][b] = spid
            self._sp_id[b][a] = spid

        self.addCovalentConstraints()
        
        self.addStructureRestraints()

        for spid, sp in self.spinpairs.items():
            a, b = sp.getAtoms()
            a, b = self._protons_num[a.getId()], self._protons_num[b.getId()]
            
            cov_score = self._get_covalent_score(a, b)
            self._sp_cov_scores[a][b] = cov_score
            self._sp_cov_scores[b][a] = cov_score

        return 1

    def setDefaultNetworkScores(self, s):

        for p in self.peaks:
            contribs = p.getContributions()
            n = len(contribs)
            [c.setNetworkScore(s/n) for c in contribs]


    # use additional distance restraints    
    def addDistanceRestraints(self):
        """
        Distance contraints
        """
        # get list of DistanceRestraints valid for NA
        restraints = []
        restraint_list =  self.iteration.getDistanceRestraints()
        
        for l, r in restraint_list.items():
            if l.getListSource()['add_to_network'] == YES:
                restraints += r
        
        if not restraints:
            return
        
        from sets import Set
        
        for r in restraints:
            for c in r.getContributions():
                for sp in c.getSpinPairs():
                    sid = sp.getId() + 1
                    self._c_id.setdefault(sid, Set())
                    self._c_id[sid].add(c)

                    if not self.spinpairs.has_key(sid):
                        self.spinpairs[sid] = sp

    def addStructureRestraints(self):

        check = {}
        vmax = self.getSettings()['v_max']
        
        for c in self.molecule.get_chains():

            residues = c.getResidues()

            atoms = [a for r in residues for a in r.getAtoms() if a.isProton() and a.getName() in ['HA', 'H']]

            for i in range(0, len(atoms)-1):
                for j in range(i+1, len(atoms)):
                    a, b = atoms[i], atoms[j]

                    id = (min(a.getId(),b.getId()), max(a.getId(),b.getId()))

                    if not check.has_key(id):
                        check[id] = 1

                        res1 =int(a.getResidue().getNumber())
                        str1 = a.getResidue().getStructure()
                        t1 = a.getName()

                        res2 = int(b.getResidue().getNumber())
                        str2 = b.getResidue().getStructure()
                        t2 = b.getName()

                        if str1 == "" or str2 == "":
                            continue

                        sep = abs(res1 - res2)
                        if sep > 4:
                            continue

                        both_H = str1 == str2 and str1[0] == 'H'
                        both_B = str1 == str2 and str1[0] == 'B'
                        if not both_B or not both_H:
                            continue

                        HA_HN = (t1 == 'HA' and t2 == 'H') or \
                                (t1 == 'H' and t2 == 'HA')

                        HN_HN = (t1 == t2) and (t1 == 'H')

                        # check if valid constraints in SS
                        
                        d = 0
                        # Sheets, dHA,HN(i,i+1)
                        if both_B and HA_HN and sep == 1:
                            d = 1

                        if both_H:
                            if HA_HN and sep <= 4:
                                d = 1

                            if HN_HN and sep <= 2:
                                d = 1

                        if d:                                
                            ##cc = CovalentConstraint(n, a, b, d)

                            a, b = self._protons_num[a.getId()], self._protons_num[b.getId()]
                            if self._sp_id[a][b] == 0:
                                self._sp_id[a][b] = -1
                            if self._sp_id[b][a] == 0:
                                self._sp_id[b][a] = -1

                            self._sp[a][b] = 2
                            self._sp[b][a] = 2

                            self._sp_cov_scores[a][b] = vmax
                            self._sp_cov_scores[b][a] = vmax

                            n+= 1

    
    def addCovalentConstraints(self):
        """
        Covalent contraints
        """
        
        dmax = self.getSettings()['distance_max']
        vmax = self.getSettings()['v_max']

        from aria.CovalentDistances import CovalentDistances
        cd = CovalentDistances()
        
        check = {}
        n = 0
        
        for c in self.molecule.get_chains():

            residues = c.getResidues()

            for r in range(len(residues)-1):

                atoms = residues[r].getAtoms() + residues[r+1].getAtoms()

                # NOESY
                atoms = [a for a in atoms if a.isProton()]

                for i in range(0, len(atoms)-1):
                    for j in range(i+1, len(atoms)):
                        aa, bb = atoms[i], atoms[j]

                        id = (min(aa.getId(),bb.getId()), max(aa.getId(),bb.getId()))
                            
                        if not check.has_key(id):
                            check[id] = 1

                            d = cd.areConnected(aa, bb)
                            
                            if d:                                
                                cc = CovalentConstraint(n, aa, bb, d)

                                a, b = self._protons_num[aa.getId()], self._protons_num[bb.getId()]
                                if self._sp_id[a][b] == 0:
                                    self._sp_id[a][b] = -1
                                if self._sp_id[b][a] == 0:
                                    self._sp_id[b][a] = -1

                                self._sp[a][b] = 2
                                self._sp[b][a] = 2
                                
                                self._sp_cov_scores[a][b] = vmax
                                self._sp_cov_scores[b][a] = vmax
                                
                                # valid also for hetero atom
                                if self._is_noesy_only:
                                    continue
                                
                                ah, bh =  aa.getHeteroAtom(), bb.getHeteroAtom()
                                if ah and bh and (ah.getType() in ['N','C'] and bh.getType() in ['N','C']) :

                                    ai, bi =  self._protons_num[ah.getId()], self._protons_num[bh.getId()]
                                    if self._sp_id[ai][bi] == 0:
                                        self._sp_id[ai][bi] = -1
                                    if self._sp_id[bi][ai] == 0:
                                        self._sp_id[bi][ai] = -1

                                    self._sp[ai][bi] = 2
                                    self._sp[bi][ai] = 2
                                
                                    self._sp_cov_scores[ai][bi] = vmax
                                    self._sp_cov_scores[bi][ai] = vmax
                                    

                                n+= 1
                                
##         # cov_score
##         for spid, sp in self.spinpairs.items():
##             a, b = sp.getAtoms()
##             d = cd.areConnected(a, b)
##             if d:
##                 map(lambda c: (c.setCovalentScore(1.)),  self._c_id[spid])
                
                
    def create_network(self):
        """
        create the network itself
        dictionnary : key = spid
                      value = Set of gammas
        """


        if self.anchoring is not None:
            return
        
        self.message('Creating network ...')
        
        from sets import Set
        
        self.anchoring = {}

        #t1 = clock()

        for spid, sp in self.spinpairs.items():


            a, b = sp.getAtoms()
            sa, sb = a.getSegid(), b.getSegid()
            a, b = self._protons_num[a.getId()], self._protons_num[b.getId()]

            # dim0
            r = self._residues_num[sa][a]
            res_bound = []
            for i in range(r-1, r+2):
                if self._residues_id[sa].has_key(i):
                    res_bound += self._residues_id[sa][i]

            x = numpy.take(self._sp, res_bound, axis = 0)

            both_0 =  x[:,a] * x[:,b]
            x0 = [res_bound[i] for i in numpy.flatnonzero(both_0)]


            # dim1
            r = self._residues_num[sb][b]
            res_bound = []
            for i in range(r-1, r+2):
                if self._residues_id[sb].has_key(i):
                    res_bound += self._residues_id[sb][i]

            x = numpy.take(self._sp, res_bound, axis = 1)

            both_1 =  x[a,:] * x[b,:]
            x1 =  [res_bound[i] for i in numpy.flatnonzero(both_1)]

            x12 = Set(x0).union(x1)
            
            self.anchoring[spid] = x12


        self.message("Done.")

                                                                                                            

    def _get_covalent_score(self, id_a, id_b):

        """
        score according to covalent structure (a, b, are two atoms
          
             { Vmax if covalent constraint
         S = { Vmin if intraresidual /sequential connectivity
             { 0 if long-range connectivity

        """
        # argument : contribution ? => then get max distance from contribution's spinpairs (use ISPA Model)
        # a spin pairs ?
        # 2 atoms

        vmin = self.getSettings()['v_min']
        vmax = self.getSettings()['v_max']

        if self._sp[id_a][id_b] == 2:
            covalent_score = vmax

        else:
            
            if self._isSequential(id_a, id_b):
                covalent_score = vmin

            else:
                covalent_score = 0.

        return  covalent_score
                

    def _heaviside(self, x):

        if x < 0:
            return 0.
        elif x == 0:
            return .5
        else:
            return 1.

    def _isSequential(self, id_a, id_b):
        sa, sb = self._protons_id[id_a].getSegid(), self._protons_id[id_b].getSegid()
        if sa <> sb:
            return 0
        else:
            return abs(self._residues_num[sa][id_a] - self._residues_num[sb][id_b]) <= 1


    def _sumContribScore(self):

        self._sp_sum_scores = {}
        
        for spid, contribs in self._c_id.items():
            s = numpy.sum([c.getScore()/len(c.getSpinPairs()) for c in contribs])
            self._sp_sum_scores[spid] = s
                
            

    def updateContributionsNetworkScores(self):

        """
        calulate network_score for each contribution and update network_score
        """

        contribs_scores = {}

        #t = clock()
        self._sumContribScore()

        #t = clock()

        v_min = self.getSettings()['v_min']
        
        for k, gammas in self.anchoring.items():

            sp = self.spinpairs[k]
            
            score = 0.

            a, b = sp.getAtoms()
            id_a = self._protons_num[a.getId()]
            id_b = self._protons_num[b.getId()]

            gammas = list(gammas)

            # a-g
            #g_scores_a = numpy.take(self._sp_sum_scores, numpy.take( self._sp_id[id_a,:], gammas))
            g_scores_a = [self._sp_sum_scores[x] for x in numpy.take( self._sp_id[id_a,:], gammas)]
            cov_scores_a = numpy.take(self._sp_cov_scores[id_a,:], gammas)
            nus_a = numpy.where(numpy.greater(g_scores_a, cov_scores_a), g_scores_a, cov_scores_a)
            nus_a *= numpy.greater(nus_a - v_min, 0)


            # b-g
            #g_scores_b = numpy.take(self._sp_sum_scores, numpy.take( self._sp_id[id_b,:], gammas))
            g_scores_b = [self._sp_sum_scores[x] for x in numpy.take( self._sp_id[id_b,:], gammas)]
            cov_scores_b = numpy.take(self._sp_cov_scores[id_b,:], gammas)
            nus_b = numpy.where(numpy.greater(g_scores_b, cov_scores_b), g_scores_b, cov_scores_b)
            nus_b *= numpy.greater(nus_b - v_min, 0)

            score = numpy.sum(numpy.sqrt(nus_a * nus_b))
            
            contribs = self._c_id[k]
            
            for c in contribs:
                contribs_scores.setdefault(c, [])
                contribs_scores[c].append(score)


        for c, ss in contribs_scores.items():
            c.setNetworkScore(numpy.sum(ss)/len(ss))#/len(ss)
            
        for p in self.peaks:

            contribs = p.getContributions() 

            scores = numpy.array([c.getNetworkScore() for c in contribs])
            #covalent = numpy.array([c.getCovalentScore() for c in contribs])
            #covalent = numpy.greater(covalent, 1.)
            #zero_scores_covalent = numpy.equal(scores, 0) * covalent
            
            #scores = numpy.where(zero_scores_covalent, 1., scores)
            
            sum_scores = numpy.sum(scores)
            if sum_scores > 0.:
                scores /= sum_scores

            map(lambda c,s : (c.setNetworkScore(s)), contribs, scores)
            
        #self.message("Done %5.3f" % (clock() -t))
                                    


    def updateContributionsScores(self):

        """
        calulate score of ecah contribution and update score
        """

        for p in self.peaks:

            contribs = p.getContributions() 
            #mask = [c.isInter() for c in contribs]
            
            scores = numpy.array([c.getNetworkScore() *  c.getWeight() for c in contribs])
            #numpy.putmask(scores, mask, scores * 1.5)
            
            sum_scores = numpy.sum(scores)
            if sum_scores > 0.:
                scores /= sum_scores

            map(lambda c,s : (c.setScore(s)), contribs, scores)
            
    
        #self.message("Done %5.3f" % (clock() -t))

    def dump_text(self):
        
        settings = None
        peak_list = self.peaks
        
        itn = self.iteration.getNumber()
        infra = self.project.getInfrastructure()
        
        import os
        from aria.Protocol import REPORT_NOE_RESTRAINTS
        
        path = infra.get_iteration_path(itn)        
        filename = os.path.join(path, REPORT_NOE_RESTRAINTS + '.network')

        pickler = NetworkAnchoringTextPickler(settings)
        pickler.dump_network(peak_list, filename, gzip = 0)
        self.message('Network-Anchoring scores (text) written (%s).' % filename)

    def dump_ps(self):
        
        itn = self.iteration.getNumber()
        infra = self.project.getInfrastructure()
        
        import os
        from aria.Protocol import REPORT_NOE_RESTRAINTS
        
        path = infra.get_iteration_path(itn)
        path = os.path.join(path, 'graphics/network')
        
        np = NetworkPsPickler(self)
        try:
            np.plot(path)
        except Exception, msg:
            import aria.tools as tools
            self.warning(tools.last_traceback())
            msg = 'Error during creation of %s.network.' % REPORT_NOE_RESTRAINTS
            self.warning(msg)
    

    def _dump_scores(self, old_weights):
        
        ## save scores
        s = ""
        n = 0
        for p in self.peaks:
            pnetscores = self.getPeakNetScores(p)
            for c in p.getContributions():
                s += "NETWORK : I %4d %5d OW %5.3f W %5.3f N %5.3f S %5.3f Nres %5.3f Nat %5.3f\n" \
                     %(p.getId(), c.getId(), old_weights[n], c.getWeight(), c.getNetworkScore(), \
                       c.getScore(), pnetscores['residue'], pnetscores['atom'])

                n += 1
                
        itn = self.iteration.getNumber()
        infra = self.project.getInfrastructure()
        import os
        path = os.path.join(infra.get_iteration_path(itn), "scores.dat")

        f = open(path, 'w')
        f.write(s)
        f.close()

        s = ''
        for k, v in self.residue_score.items():
            s += "%d %d %.4f\n" % (k[0],k[1], v)

        path = os.path.join(infra.get_iteration_path(itn), "res_scores.dat")

        f = open(path, 'w')
        f.write(s)
        f.close()              
            

    def getPeakNetScores(self, p):
        
        score = {'residue' : 0.,
                 'atom' : 0.}
        
        for c in p.getContributions():
            res = [0,1]
            for a in res:
                res[a] = c.getSpinSystems()[a].getAtoms()[0].getResidue()
            
            score['residue'] += self.getResNetScore(res) * c.getScore()
            score['atom'] += c.getNetworkScore()/len(c.getSpinPairs()) * c.getScore()

        return score    
                                   
    def getResNetScore(self, residues):

        residues.sort(lambda a,b: cmp(a.getNumber(), b.getNumber()))
        key = tuple(residues)
        
        #r1, r2 = residues[0].getNumber(), residues[1].getNumber()
        #key = (min((r1, r2)), max((r1, r2)))
        return self.residue_score[key]
                                                                
    
    def analyze(self):
        

        """
        Analyse contribution scores and remove non valable ones
        """

        self.message('Analyzing ...')

        self.result = {}
        result = {}
        
        # compute net score per residue pairs
        self.residue_score = {}

        for spid, sp in self.spinpairs.items():
            a, b = sp.getAtoms()
            r1, r2 = a.getResidue(), b.getResidue()
            #sa, sb = a.getSegid(), b.getSegid()
            #a, b = self._protons_num[a.getId()], self._protons_num[b.getId()]
            #r1, r2 = self._residues_num[sa][a], self._residues_num[sb][b]
            
            key = [r1, r2]
            key.sort(lambda a,b: cmp(a.getNumber(), b.getNumber()))
            #key = (min((r1, r2)), max((r1, r2)))
            key = tuple(key)
            self.residue_score.setdefault(key, 0.)
            sc = max([c.getNetworkScore() for c in self._c_id[spid]])
            self.residue_score[key] += sc
            

        contribs = [c for p in self.peaks for c in p.getContributions()]
        scores = [c.getScore() for c in contribs]
        total = len(contribs)
        #eliminated = [c for c in contribs if c.getScore() <= 0.]
        eliminated = numpy.sum(numpy.less_equal(scores, 0.))

        self.result['total'] = total
        self.result['eliminated'] = eliminated
        self.result['ratio'] = self.result['eliminated']*100./float(total)

        
        ## SET SCORE as Weight

        old_weights = [c.getWeight() for c in contribs]
        ## save scores
        #self._dump_scores(old_weights)


        
        [c.setWeight(c.getScore()) for c in contribs]

                
        ####################################################
        
        # FILTER PEAKS according to Nres et Natom
        # First rule : <Nres>p >= Nhigh
        # OR
        # second rule : <Nres>p >= Nres_min AND <Natom>p >= Natom_min

        #Nhigh = 4.
        #Nres_min = 1.
        #Nat_min = 0.25

        s = self.getSettings()
        
        for p in self.peaks:

            res_score = self.getPeakNetScores(p)
            p._network = res_score
            
            if p.getReferencePeak().isReliable():
                continue
            
##             if not p.isAmbiguous() and p.getActiveContributions() and p.getActiveContributions()[0].isInter():
##                 continue
            
            if not (res_score['residue'] >= s['high_residue_threshold'] or         
               (res_score['residue'] >=  s['min_residue_threshold'] and res_score['atom'] >= s['min_atom_threshold'])):
                p.isActive(0)
    
        
    def update_scores(self):

        self.setDefaultNetworkScores(1.)

        #[ c.setScore(c.getNetworkScore() * c.getWeight()) for p in self.peaks for c in p.getContributions()]
        
        self.updateContributionsScores()

        n = 0

        while n < 3:    

            self.message("Round %d ..." % n)
            t = clock()
            self._round = n
            self.updateContributionsNetworkScores()
            self.updateContributionsScores()
            self.debug('Time: %ss' % str(clock() - t))
            
            n += 1            
            
        
    def run(self, iteration):

        """
        run network anchoring.
        """

        self.iteration = iteration
        self.peaks = iteration.getPeakList()

        restraints = []
        restraint_list =  self.iteration.getDistanceRestraints()
        
        for l, r in restraint_list.items():
            if l.getListSource()['filter_contributions'] == YES and \
                   l.getListSource()['run_network_anchoring']  == YES :
                
                restraints += r

        self.peaks += restraints

        self._is_noesy_only = 1 
        # check if we have non H-H pairs
        for p in self.peaks:
            contributions = p.getActiveContributions()
            if not contributions:
                continue
            atom1, atom2 = contributions[0].getSpinPairs()[0].getAtoms()
            if not atom1.isProton() and not atom1.isProton():
                self._is_noesy_only = 0
                break
        
        from aria.Singleton import ProjectSingleton
        self.project = ProjectSingleton()
        self.molecule = self.project.getMolecule()
            
        # 1) initalize       
        done = self.setup()
        
        if not done:
            s = 'Aborting. No valid peaks or restraints.' 
            self.warning(s)
            return
        
        # 2') create network

        t1 = clock()
        
        self.create_network()
        
        self.debug('Time: %ss' % str(clock() - t1))
        
        # 2) assign  network scores to contributions

        self.update_scores()

        # 4)  Analysis
        t1 = clock()
        self.analyze() 
        s = 'Done. %(eliminated)d/%(total)d (%(ratio)5.2f %%) assignment possibilities removed.\n' 
        
        self.message(s % self.result)
        self.debug('Time: %ss' % str(clock() - t1))
        
        # 5) logs
        self.dump_text()

        self.dump_ps()

        #self.halt()

class NetworkXMLPickler(XMLBasePickler):

    def _xml_state(self, x):
        e = XMLElement()

        e.enabled = x['enabled']
        e.high_residue_threshold = x['high_residue_threshold']
        e.min_residue_threshold = x['min_residue_threshold']
        e.min_atom_threshold = x['min_atom_threshold']
        
        return e

    def load_from_element(self, e):
        s = NetworkSettings()

        s['enabled'] = str(e.enabled)
        s['high_residue_threshold'] = float(e.high_residue_threshold)
        s['min_residue_threshold'] = float(e.min_residue_threshold)
        s['min_atom_threshold'] = float(e.min_atom_threshold)        
        
        return s

NetworkSettings._xml_state = NetworkXMLPickler()._xml_state


## TEST
if __name__ == '__main__':

    molecule_file = '~/devel/aria2.2_release/test/run3/data/sequence/hrdc.xml'
    ariapeaks_file='~/devel/aria2.2_release/test/run3/structures/it0/noe_restraints.pickle'
    project_file = '~/devel/aria2.2_release/test/werner2.xml'

    # read molecule
    import aria.AriaXML as AriaXML
    pickler = AriaXML.AriaXMLPickler()
    molecule =  pickler.load(molecule_file)

    
    # read pickled ariapeak list
    from aria.tools import Load
    aria_peaks = Load(ariapeaks_file)
    

    project = pickler.load(project_file)
    project.ccpn_data_sources = ()
    project.read_molecule()

    ns = project.getProtocol().getSettings()['iteration_settings'][0]['network_anchoring_settings']
    N = NetworkAnchoring(ns)
    
    class it:

        def __init__(self, peaks, n):

            self.peaks = peaks
            self.n = n
            
        def getPeakList(self):
            return self.peaks

        def getNumber(self):
            return self.n

    N.run(it(aria_peaks, 0))

    #N.dump()
    

    
