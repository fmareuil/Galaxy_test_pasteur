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
$Date: 2010/03/23 15:27:16 $
"""


import sys
from string import join
from distutils.version import LooseVersion

NUMPY_MIN_VERSION = LooseVersion('1.0')
PYLAB_NUMPY_MIN_VERSION = LooseVersion('0.87.7')
SCIPY_NUMPY_MIN_VERSION = LooseVersion('2.7.2')
CCPN_MIN_VERSION = LooseVersion('2.0')

def check_python():

    print 'Python version   ',
    
    version = float(sys.version[:3])
    
    if version < 2.4:
        print 'Python version 2.4 or higher required.'
        print 'Current version is', sys.version[:5]
    else:
        print 'ok.'

def check_numeric():

    print 'Numpy module   ', 

    msg = ""
    failed = 0
    
    try:
       import numpy
       version = LooseVersion(numpy.__version__)

       if version >= NUMPY_MIN_VERSION:
           msg =  'ok (numpy)'
           print msg
           return 'numpy'
       
       else:
           msg += 'version > %s required for Numpy' % str(NUMPY_MIN_VERSION)
           
    except:        
        msg += 'could not import Numpy module.'

    print msg
    return None

def check_numeric_slice():

    # BARDIAUX
    msg_NUMERIC = "\nThe version of Numeric (%s) is known to be incompatible with ARIA.\nConsider reverting to a more stable version (like 23.8).\n"
    msg_NUMPY   = "\nThe version of numpy (%s) is known to be incompatible with ARIA.\nConsider reverting to a more stable version.\n"
    
    msg = {'numpy'   : msg_NUMPY,
           'numeric' : msg_NUMERIC}
    
    numerix = 'numpy'
    
    try:
        from numpy import ones
        from numpy import __version__ as NUMERIC_VERSION        
    except:
	print "No numpy found."
    
    if not len(ones(10)[2:]):
        print msg[numerix] % NUMERIC_VERSION

## def check_numeric_slice():

##     msg = "\nWARNING: This version of Numeric (%s) is known to be incompatible with ARIA.\nConsider reverting to a more stable version (like 23.8).\n"
    
##     from Numeric import ones, __version__ as NUMERIC_VERSION
##     if not ones(10)[2:]:
##         print msg % NUMERIC_VERSION
    
    
def check_tix():

    print 'Tkinter and Tix modules...'

    failed_modules = []

    try:
        import Tkinter
        print 'Tkinter imported (Tk version %.3f)' % Tkinter.TkVersion
    except:
        failed_modules.append('Tkinter')
        
    try:
        import _tkinter
    except:
        failed_modules.append('_tkinter')
        
    try:
        import Tix
        Tix.Tk()
        print 'Tix imported.'
    except:
        failed_modules.append('Tix')

    if len(failed_modules):
        print 'could not import module(s) ' + join(failed_modules, '/')

def check_scientific():

    print 'ScientificPython module   ', 
    result = 0
    try:
        import Scientific.IO.PDB
        print 'ok.'
        result = 1
    except:
        print 'could not import ScientificPython module.'

    return result

def check_pylab():

    print 'Matplotlib module (optional)  ', 
    result = 0
    try:
        import matplotlib.pylab
        print 'ok.'
        result = 1
    except:
        print 'could not import Matplotlib module.'

    return result

def check_ccpn():

    print '\nCCPN distribution:',

    missing = []

    try:
        import ccpnmr
        from memops.general.Constants import currentModelVersion
        ccpn_version = LooseVersion(str(currentModelVersion))
        if ccpn_version < CCPN_MIN_VERSION:
            print ' failed.'
            print 'CCPN version >= 2.0 required (current version %s).' % ccpn_version
            return 
        else:
            print 'ok.'
    except:
        missing.append('ccpnmr')
        print
        
    try:
        import ccpnmr.format
        print 'Format converter: ok.'
    except:
        missing.append('ccpnmr.format')

    try:
        import ccpnmr.analysis
        print 'Anaysis: ok.'
    except:
        missing.append('ccpnmr.analysis')
    
    if missing:
        print 'Could not import the following modules:',
        print ', '.join(missing)
        print 'This does not matter, as long as you do not intend to use the CCPN data model. Otherwise, please make sure that the environment variables CCPNMR_TOP_DIR and PYTHONPATH are set correctly.'


if __name__ == '__main__':

    print 'Checking Python environment...\n'

    check_python()
    numeric = check_numeric()
    check_numeric_slice()
    #scientific = check_scientific()
    scientific = 1
    check_tix()
    check_ccpn()
    pylab = check_pylab()
    

    # some infos
    if not scientific or not pylab:
        if numeric == 'numpy':
            print "\nNOTE:"
            if not scientific:
                print 'Using ScientificPython with Numpy requires version >= %s of ScientificPython' % str(SCIPY_NUMPY_MIN_VERSION)
            if not pylab:
                print 'Using Matplotlib with Numpy requires version >= %s of Matplotlib' % str(PYLAB_NUMPY_MIN_VERSION)
        
