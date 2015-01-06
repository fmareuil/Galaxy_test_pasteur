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


from aria.ariabase import *
        
class Datum(AriaBaseClass):

    def __init__(self, value, error = None):
        """
        Initialize with a value for the measurement and an error bound.
        """

        check_type(value, FLOAT, FLOAT64, NONE)
        check_type(error, FLOAT, FLOAT64, NONE)
        
        AriaBaseClass.__init__(self)
        
        self.__value = value
        self.__error = error

    def __getitem__(self, index):
        if index == 0:
            return self.getValue()
        elif index == 1:
            return self.getError()
        else:
            raise IndexError, 'Invalid index. Must be 0 or 1.'

    def getValue(self):
        return self.__value

    def getError(self):
        return self.__error

    def __len__(self):
        return 2
        
    def __str__(self):

        class_name = self.__class__.__name__
        
        return '%s(val=%s, err=%s)' % (class_name, str(self[0]), str(self[1]))

    __repr__ = __str__

class DatumXMLPickler:

    order = ['value', 'error']

    def _xml_state(self, x):

        from aria.xmlutils import XMLElement

        e = XMLElement(tag_order = self.order)

        e.value = x.getValue()
        e.error = x.getError()

        return e

    def load_from_element(self, e):

        ## if e.value is '', val is set to None
        ## otherwise, Datum will validate its
        ## arguments.

        try:
            val = eval(e.value)
        except:
            val = None

        try:
            err = eval(e.error)
        except:
            err = None

        return Datum(val, err)

class ChemicalShift(Datum):

    def __init__(self, *args, **kw):
        Datum.__init__(self, *args, **kw)

    # BARDIAUX 2.2
    def __hash__(self):
        return id(self)
        
    def __eq__(self, other):

        return self.getValue() == other.getValue() and \
               self.getError() == other.getError()

class ChemicalShiftXMLPickler(DatumXMLPickler):

    def load_from_element(self, e):

        val, err = DatumXMLPickler.load_from_element(self, e)
        return ChemicalShift(val, err)

Datum._xml_state = DatumXMLPickler()._xml_state
ChemicalShift._xml_state = ChemicalShiftXMLPickler()._xml_state

if __name__ == '__main__':

    d = Datum(5.,1.)
    cs1 = ChemicalShift(7.0, 0.5)
    cs2 = ChemicalShift(6.0)
