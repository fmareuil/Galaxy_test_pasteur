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



def Dump(this, filename, mode = 'w', as_string = 0, gzip = 0):
    """
    Dump(this, filename, gzip = 0)
    Supports also '~' or '~user'.
    """

    import os, cPickle

    ## otherwise we might get problems
    ## when pickling numeric arrays

    try:
        from numpy import array
    except:
        pass
    
    if as_string:
        return cPickle.dumps(this, 1)
    
    filename = os.path.expanduser(filename)

    if not mode in ['w', 'a']:
        raise "mode has to be 'w' (write) or 'a' (append)"

    if gzip:
        open_func = gzip_open
    else:
        open_func = open
        
    f = open_func(filename, mode)
    
    cPickle.dump(this, f, 1)
        
    f.close()

def Load(filename, gzip = 0):
    import cPickle, os

    filename = os.path.expanduser(filename)

    if gzip:
        open_func = gzip_open
    else:
        open_func = open
    
    f = open_func(filename)

    objects = []

    eof = 0
    n = 0
    
    while not eof:

        try:
            this = cPickle.load(f)
            objects.append(this)
            n += 1
        except EOFError:
            eof = 1

    f.close()

    if n == 1:
        return objects[0]
    else:
        return tuple(objects)

def last_traceback():
    import sys, traceback

    return ''.join(traceback.format_exception(*sys.exc_info()))

def copy_file(src, dst):

    import os, shutil

    src = os.path.expanduser(src)
    dst = os.path.expanduser(dst)

    src = os.path.abspath(src)
    dst = os.path.abspath(dst)

    if not os.path.exists(src):
        s = 'File "%s" does not exist or cannot be accessed.'
        raise IOError, s % src

    if src == dst:
        return
        
    shutil.copy(src, dst)

def cat_files(sources, dst, separator = '\n'):

    import os

    for src in sources:
        if not os.path.exists(src):
            s = 'File "%s" does not exist or cannot be accessed.'
            raise IOError, s % src

    files = [open(f).read() for f in sources]

    file = open(dst, 'w')
    file.write(separator.join(files))
    file.close()

def touch(path):

    import os

    os.system('touch %s' % path)

def wrap_string(s, length = 80, tol = 10):
    block = make_block(s, length, tol)
    return '\n'.join(block)

def indent(lines, prefix):

    tag = ' ' * len(str(prefix))

    lines[0] = prefix + lines[0]
    lines = [lines[0]] + map(lambda s, t = tag: t + s, lines[1:])

    return '\n'.join(lines)

def gzip_open(filename, mode = 'r'):
    import gzip, os

    filename = os.path.expanduser(filename)

    if mode == 'w' and filename[-3:].lower() <> '.gz':
        filename += '.gz'

    return gzip.GzipFile(filename, mode)

def as_tuple(x):
    if type(x) not in (type([]), type(())):
        return (x,)
    else:
        return tuple(x)

def string_to_segid(s):

    if not type(s) == type(''):
        m = 'Only string can be converted to segids; argument of type ' + \
            '"%s" given.' % type(s).__name__
        raise TypeError, m

    return '%4s' % s


def make_block(s, length = 80, tol = 10):
    blocks = s.split('\n')
    l = []
    for block in blocks:
        l += _make_block(block, length, tol)

    return l

def _make_block(s, length, tol):

    l = s.split(' ')
    l = [(w,' ') for w in l]

    words = []
    for ll in l:
        g = ll[0].split('/')
        g = [w+'/' for w in g]
        g[-1] = g[-1][:-1] + ' '

        words += g
    
    l = []
    line = ''

    for i in range(len(words)):
        word = words[i]

        if len(line + word) <= length:
            line += word
            
        else:
            if length - len(line) > tol:
                m = length - len(line)
                line += word[:m]
                word = word[m:]

            if len(line) > 1 and line[0] == ' ' and \
                   line[1] <> ' ':
                line = line[1:]

            l.append(line)
            line = word

    line = line[:-1]
    if len(line) > 1 and line[0] == ' ' and \
       line[1] <> ' ':
        line = line[1:]

    l.append(line)

    return l

def check_modules(modules):

## TODO: one could also check whether the associated
## file is located in the correct path
    
    import sys

    failed = []

    for m in modules:
        try:
            exec('import ' + m)
        except Exception, msg:
            failed.append(m)

    if failed:
        
        print last_traceback()

        failed = map(lambda s: s + '.py', failed)
        
        s = 'Some modules (see below) could not be imported. Either they\n' + \
            'do not exist, are errorenous or the PYTHONPATH is not set' + \
            ' correctly.'
        print s
        print '(%s)' % ', '.join(failed)
        
        sys.exit(1)

def sortProtonDimensions(spin_pair, cross_peak):

    proton1, proton2 = spin_pair.getAtoms()

    assignments1 = [atom for assi in cross_peak.getProton1Assignments()
                    for atom in assi.getAtoms()]

    assignments2 = [atom for assi in cross_peak.getProton2Assignments()
                    for atom in assi.getAtoms()]

    dimensions1 = [proton1 in assignments1, proton1 in assignments2]
    dimensions2 = [proton2 in assignments1, proton2 in assignments2]    

    protons = [None, None]

    if dimensions1.count(1) == 1 and dimensions2.count(1) == 1:

        protons[dimensions1.index(1)] = proton1
        protons[dimensions2.index(1)] = proton2

    elif dimensions1.count(1) == 1 and dimensions2.count(1) == 2:

        protons = [proton2, proton2]
        protons[dimensions1.index(1)] = proton1

    elif dimensions1.count(1) == 2 and dimensions2.count(1) == 1:

        protons = [proton1, proton1]
        protons[dimensions2.index(1)] = proton2

    elif dimensions1.count(1) == 2 and dimensions2.count(1) == 2:

        protons = [proton1, proton2]

    if None in protons:
        raise ValueError, 'spin pair could not be assigned to dimensions'

    else:
        return tuple(protons)


## BARDIAUX 2.2
def circular_permutation(s):
    s = list(s)
    n = len(s)
    cp = []
    s.extend(s)
    for i in range(0, n):
        cp.append(s[i:i+n])

    return cp
    

