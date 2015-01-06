import sys, os, shutil, tempfile

TEMPLATE1 = \
"""
<!DOCTYPE project SYSTEM "project1.0.dtd">
<project name="%(--name)s" version="1.0" author="" date="" description="" comment="" references="" working_directory="%(--temporary_path)s" temp_root="%(--temporary_path)s" run="1" file_root="Project" cache="%(--cache)s" cleanup="%(--cleanup)s">
  <data>
    <ccpn_model filename=""/>
    <molecule file="%(--seq_xml)s" format="xml" ccpn_id="">
      <linkage_definition name="%(--linkdef)s" filename=""/>
      <parameter_definition name="%(--paramdef)s" filename=""/>
      <topology_definition name="%(--topodef)s" filename=""/>
    </molecule>
"""

SPECTRUM =  \
"""    
    <spectrum enabled="yes" use_assignments="%(--flagassign)s" trust_assigned_peaks="%(--flagtrustassign)s" structural_rules="%(--structrules)s" filter_diagonal_peaks="%(--filterdiag)s" filter_unassigned_peaks="%(--filterunasspeaks)s">
      <shifts file="%(--shift_xml)s" format="xml" ccpn_id="" default_shift_error="%(--defshifterr)s"/>
      <peaks file="%(--peaks_xml)s" format="xml" ccpn_id="" peak_size="%(--peaksize)s" freq_window_proton1="%(--fwinprot1)s" freq_window_hetero1="%(--fwinhet1)s" freq_window_proton2="%(--fwinprot2)s" freq_window_hetero2="%(--fwinhet2)s">
        <lower_bound_correction value="%(--lowerbound)s" enabled="%(--flagboundL)s"/>
        <upper_bound_correction value="%(--upperbound)s" enabled="%(--flagboundU)s"/>
      </peaks>
      <experiment_data molecule_correlation_time="%(--molcorrtime)s" spectrum_mixing_time="%(--specmixtime)s" spectrometer_frequency="%(--spectfreq)s" ambiguity_type="%(--ambigtype)s"/>
    </spectrum>
"""

JCOUPLINGS = \
"""
    <jcouplings file="%(--jcouplings_tbl)s" format="tbl" ccpn_id="" enabled="%(--flagjcouplings)s" parameter_class="%(--classjs)s"/>
"""

RDCS = \
"""
    <rdcs file="%(--rdcs_tbl)s" format="tbl" ccpn_id="" enabled="%(--flagrdcs)s" parameter_class="%(--classrdcs)s"/>
"""

UNAMBIGUOUS = \
"""
    <unambiguous_distance_restraints file="%(--unamb_tbl)s" format="tbl" ccpn_id="" enabled="%(--flagunamb)s" add_to_network="no" calibrate="no" run_network_anchoring="no" filter_contributions="no"/>
"""

AMBIGUOUS = \
"""
    <ambiguous_distance_restraints file="%(--amb_tbl)s" format="tbl" ccpn_id="" enabled="%(--flagamb)s" add_to_network="no" calibrate="no" run_network_anchoring="no" filter_contributions="no"/>
"""

HBOND = \
"""
    <hbonds file="%(--hbond_tbl)s" format="tbl" ccpn_id="" enabled="%(--flaghbond)s" data_type="%(--typehbond)s"/>
"""

DIHEDRALS = \
"""
    <dihedrals file="%(--talos_tbl)s" format="tbl" ccpn_id="" enabled="%(--flagdihed)s" data_type="%(--typedihed)s"/>
"""

SSBONDS = \
"""
    <ssbonds file="%(--ssbonds_tbl)s" format="tbl" ccpn_id="" enabled="%(--flagssbond)s"/>
"""

SSBRIDGES = \
"""
    <ssbridge residue1="%(--ssbridgeresid1)s" segid1="%(--ssbridgesegid1)s" residue2="%(--ssbridgeresid2)s" segid2="%(--ssbridgesegid2)s"/>
"""

HISPATCH= \
"""
    <hispatch residue="%(--hisresid)s" segid="%(--hissegid)s" proton="%(--typhis)s"/>
"""

CISPROPATCH = \
"""
    <cispropatch residue="%(--cisproresid)s" segid="%(--cisprosegid)s"/>
"""

ZNPATCH = \
"""
    <znpatch type="%(--zntype)s" residue_zn="%(--znresdzn)s" segid_zn="%(--znsegzn)s" residue1="%(znresd1)s" segid1="%(znseg1)s" residue2="%(znresd2)s" segid2="%(znseg2)s" residue3="%(znresd3)s" segid3="%(znseg3)s" residue4="%(znresd4)s" segid4="%(znseg4)s"/>
"""

OTHERDATA = \
"""
    <other_data file="%(--otherfile)s" format="tbl" ccpn_id="" enabled="%(--flagother)s" data_type="%(--typeother)s"/>   
"""
    
INITIAL_STRUCTURE = \
"""
    <initial_structure file="%(--initstructfile)s" format="%(--formatinitstruct)s" ccpn_id="" enabled="%(--flaginitstruct)s"/>
"""

TEMPLATE2 = \
"""
    <symmetry enabled="%(--flagsym)s" method="standard" n_monomers="%(--n_monom)s" symmetry_type="%(--typesym)s" ncs_enabled="%(--ncsflag)s" packing_enabled="%(--packingflag)s"/> 
  </data>
  <structure_generation engine="cns">
    <cns local_executable="/Bis/home/fmareuil/Travail/cns_logh_32_centos3/cns_solve-1005101553.exe" keep_output="%(--cnsoutput)s" keep_restraint_files="%(--restraintfile)s" create_psf_file="%(--crepsffile)s" generate_template="%(--gentemplate)s" nonbonded_parameters="%(--nonbondedparam)s">
      <annealing_parameters>
        <unambiguous_restraints first_iteration="%(--unambigrest1it)s" k_hot="%(--ukhot)s" k_cool1_initial="%(--ukcool1init)s" k_cool1_final="%(--ukcool1fin)s" k_cool2="%(--ukcool2)s"/>
        <ambiguous_restraints first_iteration="%(--ambigrest1it)s" k_hot="%(--akhot)s" k_cool1_initial="%(--akcool1init)s" k_cool1_final="%(--akcool1fin)s" k_cool2="%(--akcool2)s"/>
        <hbond_restraints first_iteration="%(--hbondrest1it)s" k_hot="%(--hkhot)s" k_cool1_initial="%(--hkcool1init)s" k_cool1_final="%(--hkcool1fin)s" k_cool2="%(--hkcool2)s"/>
        <dihedral_restraints k_hot="%(--dkhot)s" k_cool1="%(--dkcool1)s" k_cool2="%(--dkcool2)s"/>
"""

BASE_KARPLUS = \
"""
        <karplus_restraints parameter_class="1" a="6.98" b="-1.38" c="1.72" d="-60.0" k_hot="0.0" k_cool1="0.2" k_cool2="1.0"/>
"""
KARPLUS_RESTRAINTS = \
"""
        <karplus_restraints parameter_class="%(--classjs)s" a="%(--ajs)s" b="%(--bjs)s" c="%(--cjs)" d="%(--djs)s" k_hot="%(--khotjs)s" k_cool1="%(--kcool1js)s" k_cool2="%(--kcool2js)s"/>
"""

BASE_RDC = \
"""
        <rdc_restraints parameter_class="1" method="SANI" first_iteration="0" k_hot="0.0" k_cool1="0.2" k_cool2="1.0" r="0.4" d="8.0" border_hot_initial="0.1" border_hot_final="40.0" border_cool1_initial="40.0" border_cool1_final="40.0" border_cool2_initial="40.0" border_cool2_final="40.0" center_hot_initial="0.1" center_hot_final="0.1" center_cool1_initial="10.0" center_cool1_final="10.0" center_cool2_initial="10.0" center_cool2_final="10.0"/>
"""

RDC_RESTRAINTS = \
"""
        <rdc_restraints parameter_class="%(--classrdcs)s" method="%(--metrdcs)s" first_iteration="%(--firstitrdcs)s" k_hot="%(--khotrdcs)s" k_cool1="%(--kcool1rdcs)s" k_cool2="%(--kcool2rdcs)s" r="%(--rrdcs)s" d="%(--drdcs)s" border_hot_initial="%(--bhotinirdcs)s" border_hot_final="%(--bhotfinrdcs)s" border_cool1_initial="%(--bcool1inirdcs)s" border_cool1_final="%(--bcool1finrdcs)s" border_cool2_initial="%(--bcool2inirdcs)s" border_cool2_final="%(--bcool2finrdcs)s" center_hot_initial="%(--chotinirdcs)s" center_hot_final="%(--chotfinrdcs)s" center_cool1_initial="%(--ccool1inirdcs)s" center_cool1_final="%(--ccool1finrdcs)s" center_cool2_initial="%(--ccool2inirdcs)s" center_cool2_final="%(--ccool2finrdcs)s"/>
"""

TEMPLATE3 = \
"""
        <flat_bottom_harmonic_wall m_rswitch_hot="%(--hwallmrswitchhot)s" m_rswitch_cool1="%(--hwallmrswitchcool1)s" m_rswitch_cool2="%(--hwallmrswitchcool2)s" rswitch_hot="%(--hwallrswitchhot)s" rswitch_cool1="%(--hwallrswitchcool1)s" rswitch_cool2="%(--hwallrswitchcool2)s" m_asymptote_hot="%(--hwallmasymphot)s" m_asymptote_cool1="%(--hwallmasympcool1)s" m_asymptote_cool2="%(--hwallmasympcool2)s" asymptote_hot="%(--hwallasymphot)s" asymptote_cool1="%(--hwallasympcool1)s" asymptote_cool2="%(--hwallasympcool2)s"/>
        <symmetry_restraints k_packing_hot="%(--symrestkpackhot)s" k_packing_cool1="%(--symrestkpackcool1)s" k_packing_cool2="%(--symrestkpackcool2)s" last_iteration_packing="%(--symrestlastitpack)s" k_ncs="%(--symrestkncs)s"/>
        <logharmonic_potential enabled="%(--loghflag)s" use_auto_weight="%(--loghautow)s" weight_unambig="%(--loghwunamb)s" weight_ambig="%(--loghwamb)s" weight_hbond="%(--loghwhbond)s"/>
      </annealing_parameters>
      <md_parameters dynamics="%(--mddyna)s" random_seed="%(--seed)s" tad_temp_high="%(--tadTh)s" tad_timestep_factor="%(--tadtimefact)s" cartesian_temp_high="%(--cartTh)s" cartesian_first_iteration="%(--cart1it)s" timestep="%(--timestep)s" temp_cool1_final="%(--Tcool1final)s" temp_cool2_final="%(--Tcool2final)s" steps_high="%(--stepsH)s" steps_refine="%(--stepref)s" steps_cool1="%(--stepcool1)s" steps_cool2="%(--stepcool2)s"/>
    </cns>
    <job_manager default_command="LOCAL" job_management="Synchro">
      <host enabled="yes" submit_command="csh -f" state_command="no" output_command="no" executable="/Bis/home/fmareuil/Travail/cns_logh_32_centos3/cns_solve-1005101553.exe" n_cpu="4" use_absolute_path="yes"/>
    </job_manager>
  </structure_generation>
  <protocol floating_assignment="%(--flagprotfloatassign)s">
"""

ITERATION =  \
"""
    <iteration number="%(--it)s" n_structures="%(--nstr)s" sort_criterion="%(--sortcrit)s" n_best_structures="%(--nbeststr)s" n_kept_structures="%(--nkstr)s">
      <assignment/>
      <merging method="%(--mergmet)s"/>
      <calibration relaxation_matrix="%(--calrelaxmat)s" distance_cutoff="%(--caldistcutoff)s" estimator="%(--calestim)s" error_estimator="%(--calerrestim)s"/>
      <violation_analysis sigma_mode="%(--violsigmod)s" violation_tolerance="%(--violtol)s" violation_threshold="%(--violthres)s"/>
      <partial_assignment weight_threshold="%(--partassignWthres)s" max_contributions="%(--partassignmaxcontrib)s"/>
      <network_anchoring high_residue_threshold="%(--netanchHresthres)s" enabled="%(--flagnetanch)s" min_residue_threshold="%(--netanchMresthres)s" min_atom_threshold="%(--netanchMatomthres)s"/>
    </iteration>
"""

END = \
"""
    <water_refinement solvent="%(--watrefsolv)s" n_structures="%(--watrefnstruct)s" enabled="%(--flagwatref)s" write_solvent_molecules="%(--watrefsolvmol)s"/>
  </protocol>
  <analysis>
    <structures_analysis enabled="%(--flagstructanalyse)s"/>
    <procheck executable="procheck" enabled="no"/>
    <prosa executable="prosa" enabled="no"/>
    <whatif executable="whatif" enabled="no"/>
    <clashlist executable="" enabled="no"/>
  </analysis>
  <report>
    <ccpn export_assignments="no" export_noe_restraint_list="no" export_structures="no"/>
    <molmol enabled="%(--reportfmolmol)s"/>
    <noe_restraint_list pickle_output="%(--reportfNOEpickleout)s" text_output="%(--reportfNOEtxtout)s" xml_output="%(--reportfNOExmlout)s"/>
    <spectra write_assigned="%(--reportfspectwassign)s" write_assigned_force="%(--reportfspectwassignf)s" iteration="%(--reportfspectit)s" write_unambiguous_only="%(--reportfspectunamb)s"/>
  </report>
</project>
"""

def invalid_args():

    import sys

    print "USAGE"
    sys.exit(1)

def parse_args(args):

    import getopt, string
    options =''
    long_options = ('output_init=', 'result_output=', 'name=', 'cache=', 'cleanup=', 
    'seq_xml=', 
    'linkdef=', 
    'paramdef=', 
    'topodef=', 
    'spectrums=', 'flagassign=', 'flagtrustassign=','structrules=', 'filterdiag=','filterunasspeaks=',
    'shift_xml=', 'defshifterr=', 
    'peaks_xml=', 'peaksize=','fwinprot1=','fwinhet1=','fwinprot2=','fwinhet2=',
    'lowerbound=', 'flagboundL=', 
    'upperbound=','flagboundU=', 
    'molcorrtime=','specmixtime=','spectfreq=','ambigtype=',
    'Jrestraints=','jcouplings_tbl=','flagjcouplings=','classjs=','ajs=','bjs=','cjs=','djs=','khotjs=','kcool1js=', 'kcool2js=',
    'RDCrestraints=','rdcs_tbl=','flagrdcs=','classrdcs=','metrdcs=','firstitrdcs=','khotrdcs=','kcool1rdcs=','kcool2rdcs=','rrdcs=','drdcs=','bhotinirdcs=','bhotfinrdcs=','bcool1inirdcs=','bcool1finrdcs=','bcool2inirdcs=','bcool2finrdcs=','chotinirdcs=','chotfinrdcs=','ccool1inirdcs=','ccool1finrdcs=','ccool2inirdcs=','ccool2finrdcs=',
    'Urestraints=','unamb_tbl=','flagunamb=',
    'Arestraints=','amb_tbl=','flagamb=',
    'Hrestraints=','hbond_tbl=','flaghbond=','typehbond=',
    'Drestraints=','talos_tbl=','flagdihed=','typedihed=',
    'flagsym=','n_monom=','typesym=','ncsflag=','packingflag=',
    'SS1restraints=','ssbonds_tbl=','flagssbond=',
    'SS2restraints=','ssbridgeresid1=','ssbridgesegid1=','ssbridgeresid2=','ssbridgesegid2=',
    'HISrestraints=','hisresid=','hissegid=','typhis=',
    'CISrestraints=','cisproresid=','cisprosegid=',
    'ZNrestraints=','zntype=','znresdzn=','znsegzn=','znresd1=','znseg1=','znresd2=','znseg2=','znresd3=','znseg3=','znresd4=','znseg4=',
    'Orestraints=','otherfile=','flagother=','typeother=',
    'Istructs=','initstructfile=','formatinitstruct=','flaginitstruct=',
    'cnsoutput=','restraintfile=','crepsffile=','gentemplate=','nonbondedparam=',
    'unambigrest1it=','ukhot=','ukcool1init=','ukcool1fin=','ukcool2=',
    'ambigrest1it=','akhot=','akcool1init=','akcool1fin=','akcool2=',
    'hbondrest1it=','hkhot=','hkcool1init=','hkcool1fin=','hkcool2=',
    'dkhot=','dkcool1=','dkcool2=',
    'hwallmrswitchhot=','hwallmrswitchcool1=','hwallmrswitchcool2=','hwallrswitchhot=','hwallrswitchcool1=','hwallrswitchcool2=','hwallmasymphot=','hwallmasympcool1=','hwallmasympcool2=','hwallasymphot=','hwallasympcool1=','hwallasympcool2=',
    'symrestkpackhot=', 'symrestkpackcool1=','symrestkpackcool2=','symrestlastitpack=','symrestkncs=',
    'loghflag=','loghautow=','loghwunamb=','loghwamb=','loghwhbond=',
    'mddyna=','seed=','tadTh=','tadtimefact=','cartTh=','cart1it=','timestep=','Tcool1final=','Tcool2final=','stepsH=','stepref=','stepcool1=','stepcool2=',
    'flagprotfloatassign=',
    'iterations=', 'it=','nstr=','sortcrit=','nbeststr=','nkstr=',
    'mergmet=',
    'calrelaxmat=','caldistcutoff=','calestim=','calerrestim=',
    'violsigmod=','violtol=','violthres=',
    'partassignWthres=','partassignmaxcontrib=',
    'netanchHresthres=','flagnetanch=','netanchMresthres=','netanchMatomthres=',
    'watrefsolv=','watrefnstruct=','flagwatref=','watrefsolvmol=',
    'flagstructanalyse=',
    'reportfmolmol=',
    'reportfNOEpickleout=','reportfNOEtxtout=','reportfNOExmlout=',
    'reportfspectwassign=','reportfspectwassignf=','reportfspectit=','reportfspectunamb=',
    'temporary_path=')
    #try:
    opt_list, args = getopt.getopt(args, options, long_options)
    #except:
    #    invalid_args()
    #print args
    #print opt_list

    d = {}
    for key, value in opt_list:
        if value == '':
            value = None
            
        d[key] = value
    for l in long_options:
       tmp = "--"+ string.split(l,"=")[0]
       if not d.has_key(tmp):
           d[tmp] = ''
       if tmp != '--typesym':
         if d[tmp] == "None" or d[tmp] == None:
           d[tmp] = ''   
    return d

def writeboucle(liste, TEMPLATE):
  for i in string.split(liste, ","):
    tmp = string.split(i)
    if len(tmp) > 0:
        dicospec = parse_args(tmp)
        NEWTEMPLATE = TEMPLATE % dicospec
	Project_conv.write(NEWTEMPLATE)

####MAIN####
import string
dico = parse_args(sys.argv[1:])
Project_conv = open(dico["--output_init"],"w")
os.chdir(dico["--temporary_path"])
tmpd = tempfile.mkdtemp(prefix=os.getcwd()+"/tmp")
dico["--temporary_path"] = tmpd
os.chdir(dico["--temporary_path"])

TEMPLATE1 = TEMPLATE1 % dico
Project_conv.write(TEMPLATE1)

writeboucle(dico["--spectrums"],SPECTRUM)
comp = 1
for j in string.split(dico["--Jrestraints"], ","):
    tmp = string.split(j)
    if len(tmp) > 0:
        dicospec = parse_args(tmp)
        dicospec["--classjs"] = comp
        NEWJCOUPLINGS = JCOUPLINGS % dicospec
	Project_conv.write(NEWJCOUPLINGS)
        comp += 1

comp = 1
for j in string.split(dico["--RDCrestraints"], ","):
    tmp = string.split(j)
    if len(tmp) > 0:
        dicospec = parse_args(tmp)
        dicospec["--classrdcs"] = comp
        NEWRDCS = RDCS % dicospec
	Project_conv.write(NEWRDCS)
        comp += 1
writeboucle(dico["--Urestraints"],UNAMBIGUOUS)
writeboucle(dico["--Arestraints"],AMBIGUOUS)
writeboucle(dico["--Hrestraints"],HBOND)
writeboucle(dico["--Drestraints"],DIHEDRALS)
writeboucle(dico["--SS1restraints"],SSBONDS)
writeboucle(dico["--SS2restraints"],SSBRIDGES)
writeboucle(dico["--HISrestraints"],HISPATCH)
writeboucle(dico["--CISrestraints"],CISPROPATCH)
writeboucle(dico["--ZNrestraints"],ZNPATCH)
writeboucle(dico["--Orestraints"],OTHERDATA)
writeboucle(dico["--Istructs"],INITIAL_STRUCTURE)

TEMPLATE2 = TEMPLATE2 % dico				
Project_conv.write(TEMPLATE2)

comp = 1
for j in string.split(dico["--Jrestraints"], ","):
    tmp = string.split(j)
    if len(tmp) > 0:
        dicospec = parse_args(tmp)
        dicospec["--classjs"] = comp
        NEWKARPLUS_RESTRAINTS = KARPLUS_RESTRAINTS % dicospec
	Project_conv.write(NEWKARPLUS_RESTRAINTS)
        comp += 1

if len(dico["--Jrestraints"].split()) == 0:
    Project_conv.write(BASE_KARPLUS)

comp = 1
for j in string.split(dico["--RDCrestraints"], ","):
    tmp = string.split(j)
    if len(tmp) > 0:
        dicospec = parse_args(tmp)
        dicospec["--classrdcs"] = comp
        NEWRDC_RESTRAINTS = RDC_RESTRAINTS % dicospec
	Project_conv.write(NEWRDC_RESTRAINTS)
        comp += 1

if len(dico["--RDCrestraints"].split()) == 0:	
    Project_conv.write(BASE_RDC)
     
TEMPLATE3 = TEMPLATE3 % dico	
Project_conv.write(TEMPLATE3)

comp = 0
for j in string.split(dico["--iterations"], ","):
    tmp = string.split(j)
    if len(tmp) > 0:
        dicospec = parse_args(tmp)
        dicospec["--it"] = comp
        NEWITERATION = ITERATION % dicospec
	Project_conv.write(NEWITERATION)
        comp += 1

END = END % dico	
Project_conv.write(END)
Project_conv.close()

command = "python -W ignore::DeprecationWarning -O /Bis/home/fmareuil/Travail/Galaxy/galaxy-dist/tools/mytools/ARIA/aria2.py -s %s" % (dico["--output_init"])
os.system(command)

command = "python -W ignore::DeprecationWarning -O /Bis/home/fmareuil/Travail/Galaxy/galaxy-dist/tools/mytools/ARIA/aria2.py %s" % (dico["--output_init"])
os.system(command)
if os.path.isdir("%s/run1/structures/refine" % dico["--temporary_path"]):
  command = "tar -P -cvf %s %s/run1/structures/it%i %s/run1/structures/refine" % (dico["--result_output"], dico["--temporary_path"], comp-1,  dico["--temporary_path"])
else:
  command = "tar -P -cvf %s %s/run1/structures/it%i " % (dico["--result_output"], dico["--temporary_path"], comp-1)
os.system(command)
shutil.rmtree("%s/run1" % dico["--temporary_path"])
