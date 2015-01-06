#!/usr/bin/env python

import os, sys, string, urllib2, smtplib, stat
sys.path.insert( 0, os.path.dirname( __file__ ) )
from common import display, get, post, submit
from xml.etree import ElementTree as ET

def main():
    try:
        key = sys.argv[1]
        url = sys.argv[2]
    except IndexError:
        print >> sys.stderr, 'usage: %s key url' % os.path.basename( sys.argv[0] )
        sys.exit( 1 )  
    return tool_list( key, url)

def getcommand(dico):
   listinterdit = ['PYTHONPATH="/mount/biotools/lib/python2.6/site-packages/"','#if','##', 'stderr_wrapper.sh', '$db_opts.db_opts_selector', '==', "'available_db'", ":"]
   for key in dico:
     if "galaxy.web.pasteur.fr" in key or "toolshed.g2.bx.psu.edu" in key or "testtoolshed.g2.bx.psu.edu" in key:
        base_dir = "../../../shed_tools/%s" % os.path.dirname(os.path.dirname(key))
        basenam = os.path.basename(os.path.dirname(os.path.dirname(key)))       
        for version in os.listdir(base_dir):
            if "urgi_tools" in key:
                base_version = os.path.join(base_dir,version,basenam,"repet_pipe/SMART/galaxy")
            elif "picard_pasteur_wrapper" in key:
                base_version = os.path.join(base_dir,version,basenam,basenam)
            elif "bcftools" in key:
                base_version = os.path.join(base_dir,version,basenam,"bcftools-3182c7fac413")
            elif "samtools_idxstats" in key:
                base_version = os.path.join(base_dir,version,basenam,"tools/samtools_idxstats")
            elif "brad-chapman" in key:
                base_version = os.path.join(base_dir,version,basenam,basenam)
            else:
                base_version = os.path.join(base_dir,version,basenam)
                #if os.path.exists(base_dir):
                #tagxmlfound = 0
            for namedoc in [x.encode('utf-8') for x in os.listdir(base_version)]:
               if len(dico[key]) == 3: 
                    tool_id = None
                    if "xml" in namedoc and not "xml~" in namedoc and not namedoc.startswith(".") and not namedoc.endswith(".sample") and not namedoc == "datatypes_conf.xml" and not namedoc == "tool_dependencies.xml":
                        #tagxmlfound = 1
                        doc = ET.parse( os.path.join(base_version,namedoc))
                        root = doc.getroot()
                        #if root.tag == "tool":
                        if root.attrib["id"] == os.path.basename(os.path.dirname(key)):
                            #print root.attrib["id"]
                            for rootchild in root.getchildren():
                                if  rootchild.tag == "command" :
                                    for command in rootchild.text.split():
                                          if command not in listinterdit:
                                            interpretor = None
                                            try: 
                                                interpretor = rootchild.attrib["interpreter"]
                                                dico[key].append(command + "\t" + interpretor)
                                                break
                                            except:
                                                dico[key].append(command)
                                                break

   for key in  dico :
        if len(dico[key]) == 4:
          if "galaxy.web.pasteur.fr" in key:
            print "Pasteur tool\t%s\t%s\t%s\t%s\t%s" % (os.path.basename(os.path.dirname(key)), dico[key][0], dico[key][1], dico[key][2], dico[key][3])
          elif "toolshed.g2.bx.psu.edu" in key:
            print "Toolshed tool\t%s\t%s\t%s\t%s\t%s" % (os.path.basename(os.path.dirname(key)), dico[key][0], dico[key][1], dico[key][2], dico[key][3])
          else:
            print "Bizarre tool\t%s\t%s\t%s\t%s\t%s" % (os.path.basename(os.path.dirname(key)), dico[key][0], dico[key][1], dico[key][2], dico[key][3])
        
        else:
          print "Xml not in toolsheds tool\t%s\t%s\t%s\t%s" % (key, dico[key][0], dico[key][1], dico[key][2])
                        #       if ( rootchild.tag == "command") and goodtool == 1:
                        #           index = 0
                        #           print rootchild.text.split()[index]

def tool_list(key, url):
    try:
        tools = display(key, url, return_formatted=False)
    except Exception, e:
        print >> sys.stderr, ("Error during the display API galaxy command : %s" % e)
        sys.exit( 1 )
    
    
    not_good_command = ['stderr_wrapper.sh']
    panel = None    
    #num = 0
    dict_ids = {}
    for i in tools:
        for elem in i["elems"]:
            if elem.has_key("panel_section_name"): 
                panel_section = elem["panel_section_name"] 
            else: 
                panel_section = None
            if elem.has_key("id"):
                if "toolshed.g2.bx.psu.edu" in elem["id"] or "testtoolshed.g2.bx.psu.edu" in elem["id"]:
                    id = "Toolshed tool\t%s" % os.path.basename(os.path.dirname(elem["id"]))
                elif "galaxy.web.pasteur.fr" in elem["id"]:
                   id = "Pasteur tool\t%s" % os.path.basename(os.path.dirname(elem["id"]))
                   #display(key, os.path.join(url + "%s?io_details=True" % elem['id']))
                else:
                    id = "Galaxy tool\t%s" % elem["id"]
            else: 
                id = None
            if elem.has_key("version"): 
                version = elem["version"] 
            else: 
                version = None
            if elem.has_key("name"): 
                name = elem["name"] 
            else: 
                name = None

            if panel_section and id and version and name:
                dict_ids[elem["id"]] = [version, name, panel_section]

    return dict_ids        
            
if __name__=="__main__":
    dict_tools = main()
    getcommand(dict_tools)
