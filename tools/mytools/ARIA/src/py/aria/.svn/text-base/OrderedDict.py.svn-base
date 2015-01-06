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


## TODO: get rid of UserDict
UserDict = dict

class OrderedDict(UserDict):

    def __init__(self, order = None):

        UserDict.__init__(self)
        self.order = order

    def keys(self):

        if self.order is not None:
            return self.order
        else:
            return UserDict.keys(self)

    def values(self):
        return map(lambda k, s = self: s[k], self.keys())

    def items(self):
        return map(lambda k, s = self: (k, s[k]), self.keys())

    def __setitem__(self, key, value):
        
        if self.order is None:
            self.order = []

        if key not in self.order:
            self.order.append(key) 

        UserDict.__setitem__(self, key, value)

    def __delitem__(self, name):
        if self.order is not None:
            if name in self.order:
                self.order.remove(name)

        UserDict.__delitem__(self, name)

