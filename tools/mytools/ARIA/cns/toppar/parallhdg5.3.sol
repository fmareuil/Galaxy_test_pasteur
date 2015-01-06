remarks   PARAM19.SOL (solvent parameters)
remarks   ===========
remarks   available: TIPS3P and DMSO model

set echo=false end


{* TIP3P model *}
{* =========== *}
BOND HT   OT     450.0       0.9572    ! from TIPS3P geometry ( SHAKE w/PARAm)
BOND HT   HT       0.0       1.5139    ! from TIPS3P geometry ( SHAKE w/PARAm)
ANGLE HT   OT   HT      55.0     104.52    ! FROM TIPS3P geometry



{* original TIP3P for ALL interactions since we use OPLS for the protein *}
 NONBONDED  OT       .1521   3.1506    .1521   3.1506 ! ALLOW WAT 
						      !TIP3P OXYGEN PARAMETERS, adm jr., NBFIX obsolete 
 NONBONDED  HT       .0460    .4000    .0460    .4000 ! ALLOW WAT 
						      !TIP3P HYDROGEN PARAMETERS, adm jr., NBFIX obsolete 
{*            *} ! DMSO parameters Alexandre Bonvin Utrecht University
{* DMSO model *} ! GROMOS parameters from Liu et al JACS 117:4363 (1995)
{* ========== *}
BOND CDMS SDMS  300.0       1.950 
BOND SDMS ODMS  350.0       1.530

ANGLE CDMS SDMS CDMS   110.0   97.40
ANGLE CDMS SDMS ODMS   110.0  106.75

NONBONDED CDMS  0.2940  3.66  0.2940  3.66
NONBONDED SDMS  0.2384  3.56  0.2384  3.56
NONBONDED ODMS  0.0715  2.63  0.0715  2.63


set echo=true end
