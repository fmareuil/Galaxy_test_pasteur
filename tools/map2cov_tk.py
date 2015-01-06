#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Program: MAP2COV - Convert mapping data into coverage format readable by the COV2HTML interface

    Reference format supported:
        GenBank, EMBL, FASTA/GFF3

    Alignment format supported:
        SAM/BAM, ELAND, WIG

    Data type supported:
        DNA, RNA, TSS, CHip
        
    Optional Parameters
        Anonymize, Read-length, Strand-specific, Coverage Style
"""

__author__ = 'Marc Monot'
__version__= '4.4'
__email__ = 'marc.monot@pasteur.fr'


#  This file is a part of COV2HTML software
#  COV2HTML <http://sourceforge.net/projects/cov2html/> 
#  A visualization and analysis tool of Bacterial NGS data for biologists.
# 
#  ---------------------------------------------------------------------- 
#  Copyright (C) 2012 Marc Monot
#  
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 3 of the License, or
#   (at your option) any later version. 
#   <http://www.gnu.org/licenses/gpl-3.0.html>
# 
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#  ---------------------------------------------------------------------- 
#  
#  @author Marc Monot <marc.monot@pasteur.fr>
#  @author Mickael Orgeur


### PYTHON Module
# Base
import sys
import re

# Project
from _modules.functions import *
from _modules.tk import *

# Add from Python Package Index
from _modules.ordereddict import *


### INTERFACE:
if __name__ == '__main__':
    win = LaunchInterface(None)
    win.title('MAP2COV - Convert mapping data into coverage format readable by the COV2HTML interface')
    win.mainloop()

