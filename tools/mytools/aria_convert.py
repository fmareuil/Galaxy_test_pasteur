import sys, os

TEMPLATE = \
"""
<!DOCTYPE conversion SYSTEM "conversion1.0.dtd">

<!--===================== Conversion ==================================

  Template XML file for converting NMR data from various formats into
  ARIA XML format. Fill in the fields between the quotation marks "".
  Leave optional fields unchanged.

  ARIAs conversion routines supports the following formats:

  Ansig, NMRView, XEasy

  If data conversion is performed with the CCPNMR software suite (in this case
  use the ARIA command line option "convert_ccpn" instead of "convert" to
  launch the conversion) the following formats are supported:

  Ansig, NmrView, Pronto, Sparky, XEasy

-->

<conversion>

<!--==================== Project definition ===========================

  If you want ARIA to automatically create a project XML file, specify the
  desired filename and the name of your project in the fields "filename" and
  "name", respectively. During the conversion process, your data files
  are converted into ARIA XML format (the original files still exist after
  conversion) and referenced by the project file. Use ARIAs graphical user
  interface (GUI), an XML editor, or a text editor to edit the project file.

  Alternatively, you may also leave the fields empty and create an empty
  project XML file by yourself by invoking ARIA from the command line. 

  Optional fields:

  - filename

-->  

  <project name="%(--name)s">
     <output
       filename=""/>
  </project>       

<!--=================== Molecule definition ===========================

  Supported input formats are "seq" (a sequence of three letter codes) or
  "pdb". If your sequence is in PDB format, a naming convention for the
  residues and atoms needs to be specified; ARIA supports "iupac", "dyana",
  and "cns" format. Furthermore, the "molecule_segid" entry will be ignored
  and instead read from the PDB file. The type of the molecule can
  be set to "PROTEIN", "DNA" or "RNA". The additional attribute
  "first_residue_number" specifies where the residue numbering starts
  (in case of SEQ format). Fill in the attribute "molecule_name" only if you
  want to use CCPNMR conversion.

  Optional fields:

  - molecule_type        (in case of PDB input files)
  - molecule_name        
  - molecule_segid       (in case of PDB input files)
  - first_residue_number (in case of PDB input files)
  
  Conversion via CCPNMR:

  Supported formats are "Ansig", "Fasta", "NmrStar", "NmrView", 
                        "Pdb", "Sparky", "XEasy"

  When working on symmetric multimeric protein, you can specify the segids of each
  monomer by seprating the segid by an "/". Example molecule_segid="A/B"

-->

  <molecule
  
    molecule_type="%(--type)s"
    molecule_name="%(--name)s"
    molecule_segid="%(--segid)s"
    first_residue_number="%(--first_resid)s">
    
    <input
      filename="%(--seq_input)s"
      format="%(--format_seq)s"
      naming_convention=""/>
      
    <output
      filename="%(--seq_output)s"/>
  </molecule>

<!--=================== Spectrum definition ===========================

  Each spectrum is defined by a list of chemical shifts,
  and a list of NOESY cross peaks. For converting the cross peaks, you need
  to specify the dimensions that correspond to the resonances. The entry for
  the 1st and the 2nd proton dimension must be a number (1, 2, 3, or 4).
  Specify the dimension of the linked heavy nuclei in "hetero1" and "hetero2",
  respectively. If you want to specify the segment(s) for which you have
  chemical shifts or cross-peaks assignments, use the "segid"-field in the
  header of the spectrum-block. In case of multiple segments the corresponding
  segids are separated by a slash (e.g., segids="A" or segids="A/B").

  Optional fields:

  - spectrum_name
  - segids

  Data conversion via CCPNMR

  If you want to use the CCPNMR software suite and the format converter for
  data conversion, you need to specify the type of your NOESY experiment.
  The following experiments are supported:

  spectrum_type:  2D homonuclear noesy: noesy.hh
                  3D noesy (carbon):    noesy_hsqc_HCH.hhc
                  3D noesy (nitrogen):  noesy_hsqc_HNH.hhn
                  4D noesy (carbon):    noesy_hsqc_HCCH.hhcc
                  4D noesy (nitrogen):  noesy_hsqc_HNNH.hhnn
                  4D noesy (carb/nitr): noesy_hsqc_HCNH.hhcn
                  4D noesy (nitr/carb): noesy_hsqc_HNCH.hhnc
                  
  Otherwise, leave the field empty. Also, the format of chemical shifts and
  cross-peaks must be identical.

  Symmetric dimer:
  
  When working on symmetric dimers, you can specify the ambiguity level of
  your NOESY experiment. For example, if your NOESY spectrum comes from an
  asymmetric labeling experiment, you must specify spectrum_ambiguity as
  "inter". ARIA will then consider all the noes from this spectrum as intre-
  molecular ones.
  This field is not compatible with the conversion via CCPNMR.

  If the spectrum_ambiguity is set ton "inter" or "all" you must sepecify
  the segids of the involved monomers.(e.g., segids="A/B").
  
  If you are working on a monomeric protein, you can just leave this field empty
  (spectrum_ambiguity set to "intra" by default)
  
  spectrum_ambiguity:  intra    (noes involving atoms from one monomer only)
                       inter    (noes involving atoms atoms from different monomers)
                       all      (no known information, all noes are ambigous in terms of monomer)

                       


  
  Supported formats: "Ansig", "NmrView", "Pronto", "Sparky", "XEasy"

  For XEasy you may specify your cross peak assignments in a separate
  ".assign" file. In that case specify the filename (attribute "filename" in
  element "assignments"). Set filename to "" for all other cases.

-->

<!--      If you want to convert more than one spectrum, copy the block
          <spectrum> ... </spectrum>, i.e. for n spectra
          would you need n blocks. -->
"""
SPECTRUM= \
"""
  <spectrum
  
    spectrum_name="%(--spectrum_name)s"
    spectrum_type="%(--spectrum_typ)s"
    spectrum_ambiguity="%(--ambig_typ)s"
    segids="">
    
    <chemical_shifts>
      <input
        filename="%(--chemical_shift_input)s"
        format="%(--format_chemical_shift)s"/>
      <output
        filename="%(--output_chemical_shift)s"/>
    </chemical_shifts>
    
    <cross_peaks>
      <input
        filename="%(--file_cross_peaks)s"
        format="%(--format_cross_peaks)s"
        proton1="%(--num_P1)s"
        hetero1="%(--num_H1)s"
        proton2="%(--num_P2)s"
        hetero2="%(--num_H2)s"/>
        
      <output
        filename="%(--output_cross_peaks)s"/>
        
      <assignments
        filename="%(--assign_input)s"/>

      <peaks_ambiguity
        filename="%(--assign_peak_ambig)s"/>        
    </cross_peaks>
  </spectrum>
"""  
FINAL= \
"""
</conversion>
"""

def invalid_args():

    import sys

    print "USAGE"
    sys.exit(1)

def parsing(args, options):
    import string
    for h in range(len(string.split(string.join(args),"="))):
        print h, h+1


def parse_args(args):

    import getopt, string
    options =''
    long_options = ('output_convertor=', 'name=', 'type=', 'segid=', 'first_resid=', 'seq_input=', 'format_seq=', 'seq_output=', 'spectrums=', 'spectrum_name=', 'spectrum_typ=', 'ambig_typ=','chemical_shift_input=','format_chemical_shift=','file_cross_peaks=','format_cross_peaks=','num_P1=','num_H1=','num_P2=','num_H2=','assign_input=','assign_peak_ambig=','output_cross_peaks1=','output_cross_peaks2=','output_cross_peaks3=','output_cross_peaks4=','output_cross_peaks5=','output_chemical_shift1=','output_chemical_shift2=','output_chemical_shift3=','output_chemical_shift4=','output_chemical_shift5=')
    #long_options = ('output_convertor=', 'name=', 'type=', 'segid=', 'first_resid=', 'seq_input=', 'format_seq=', 'seq_output=', 'spectrums=')
    #print args
    try:
        opt_list, args = getopt.getopt(args, options, long_options)
    except:
        invalid_args()
    #print args[-1:]
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
       if d[tmp] == "None" or d[tmp] == None:
           d[tmp] = ''   
    return d
###tentative ratee avec un nombre aleatoire de fichier output
def check_temp_outdir(dic, dicspec,compt):
    dicspec["--output_cross_peaks1"] = "%s/primary_%s_output_cross_peaks%i_visible_xml" % (dic["--temporary_outdir"],dic["--idpeaks"],compt)
    dicspec["--output_chemical_shift1"] = "%s/primary_%s_output_chemical_shift%i_visible_xml" % (dic["--temporary_outdir"],dic["--idshift"],compt)
    #dicspec["--output_cross_peaks1"] = "output_cross_peaks%i" % (compt)
    #dicspec["--output_chemical_shift1"] = "output_chemical_shift%i" % (compt)
    return dicspec   
    
####MAIN####
import string
dico = parse_args(sys.argv[1:])
TEMPLATE = TEMPLATE % dico
Project_conv = open(dico["--output_convertor"],"w")
Project_conv.write(TEMPLATE)
comp = 0
for i in string.split(dico["--spectrums"], ","):
    comp += 1
    tmp = string.split(i)
    if len(tmp) > 0:
        dicospec = parse_args(tmp)
	cross_peaks = "--output_cross_peaks%i" % comp
	chemical_shift = "--output_chemical_shift%i" % comp
        dicospec["--output_cross_peaks"] = dico[cross_peaks]
        dicospec["--output_chemical_shift"] = dico[chemical_shift]
        NEWSPECTRUM = SPECTRUM % dicospec
        Project_conv.write(NEWSPECTRUM)

Project_conv.write(FINAL)
Project_conv.close()
command = "python /Bis/home/fmareuil/Travail/Galaxy/galaxy-dist/tools/mytools/ARIA/aria2.py --convert %s" % (dico["--output_convertor"])
#print command
os.system(command)


