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


from aria.ariabase import AriaBaseClass

COVALENT_DATA = 'covalent_distances.xml'



class CovalentDistances(AriaBaseClass):


    def __init__(self):

        from os.path import join

        AriaBaseClass.__init__(self)
        
        path = join(AriaBaseClass.data_path, COVALENT_DATA)
        self.load_from_xml(path)


    def load_from_xml(self, xml_file):

        import aria.xmlutils as xmlutils

        content_handler = xmlutils.XMLContentHandler()
        pickler = xmlutils.XMLPickler(content_handler)
        
        covalent_distances = pickler.load(xml_file).covalent_distances

        connections = {}
        
        for connection in covalent_distances.connection:

            residues1 = str(connection.atoms1.residues).split(',')
            residues2 = str(connection.atoms2.residues).split(',')
            atoms1 = str(connection.atoms1.names).split(',')
            atoms2 = str(connection.atoms2.names).split(',')
            d_max = float(connection.distance_max)
            separation = int(connection.separation)


            for res1 in residues1:
                for a in atoms1:
                    key1 = (res1, a)
                    connections.setdefault(key1, {})
                    for res2 in residues2:
                        for b in atoms2:
                            key2 = (res2, b)
                            connections[key1].setdefault(key2, {})
                            connections[key1][key2].update({separation : d_max})


            if separation == 0:
                for res2 in residues2:
                    for b in atoms2:
                        key1 = (res2, b)
                        connections.setdefault(key1, {})
                        for res1 in residues1:
                            for a in atoms1:
                                key2 = (res1, a)
                                connections[key1].setdefault(key2, {})
                                connections[key1][key2].update({separation : d_max})
                
        self.connections = connections

    def get_connections(self):

        return self.connections

    def areConnected(self, a, b):

        atoms = [a, b]
        atoms.sort(lambda a,b:cmp(a.getResidue().getNumber(), b.getResidue().getNumber()))
        at1, at2  = map(lambda e: e.getName(), atoms)

        if a.getSegid() <> b.getSegid():
            return None
        
        res1, res2 = map(lambda e: e.getResidue(), atoms)    
        sep = res2.getNumber() - res1.getNumber()

        if self.connections.has_key((res1.getType(), at1)):
            key1 = (res1.getType(), at1)
        else:
            if self.connections.has_key(('all', at1)):
                key1 = ('all', at1)
            else:
                return None

        d = self.connections[key1]
        
        if d.has_key((res2.getType(), at2)):
            key2 = (res2.getType(), at2)
        else:
            if d.has_key(('all', at2)):
                key2 = ('all', at2)
            else:
                return None

        d = d[key2]

        if d.has_key(sep):
            return d[sep]
        else:
            return None              
                    
                

    def error(self, exception = None, error = '', msg = None, raise_error = 1):
        if raise_error <> 1: return

        AriaBaseClass.error(self, exception, error, msg)
