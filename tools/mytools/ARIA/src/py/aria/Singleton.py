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


class Singleton:
    
    def __init__(self):

        has_single = hasattr(self.__class__, '_single')

        if has_single:
            for base in self.__class__.__bases__:
                if hasattr(base, '_single'):
                    if id(base._single) == id(self.__class__._single):
                        has_single = 0
                        break
                    
        if has_single:
            raise self.__class__._single
        else:
            self.__class__._single = self

    def __getstate__(self):
        return (self.__dict__, self.__class__._single)

    def __setstate__(self, s):
        self.__dict__ = s[0]
        self.__class__._single = s[1]
            
def instantiate(constructor, *args, **kw):
    """
    always return an instance of Singleton
    if the instance already exists, returns this instance
    otherwise a new instance of Singleton

    it uses exception handling to get the Singleton instance
    """
    
    try:
        single = constructor(*args, **kw)

    ## catches the Singleton exception and receives the instance!
        
    except constructor, instance:
        single = instance
        
    return single

def SpinPairFactory():

    from aria.Factory import SpinPairFactory as _A

    return instantiate(_A)

def SpinPairListFactory():

    from aria.Factory import SpinPairListFactory as _B

    return instantiate(_B)

def ProjectSingleton(*a, **b):

    from aria.Project import ProjectSingleton as X

    if b.get('__new_instance__', 0) and hasattr(X, '_single'):
        del X._single

    if b.has_key('__new_instance__'):
        del b['__new_instance__']

    return instantiate(X, *a, **b)

def AtomFactory(*a, **kw):

    from aria.Factory import AtomFactory as f

    if kw.get('__new_instance__', 0) and hasattr(f, '_single'):
        del f._single
        
    return instantiate(f, *a)
 
