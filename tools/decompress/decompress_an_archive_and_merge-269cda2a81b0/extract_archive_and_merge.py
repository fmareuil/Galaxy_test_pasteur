#!/usr/bin/env python
"""


Created by Cyril MONJEAUD
Cyril.Monjeaud@irisa.fr

And with the help of Anthony Bretaudeau for some stuff with bz2.

"""


import argparse, os, sys, subprocess, tempfile, shutil, gzip, zipfile, tarfile, gzip, bz2
import glob 
from galaxy import eggs
from galaxy import util
from galaxy.datatypes.checkers import *

def stop_err( msg ):
    sys.stderr.write( '%s\n' % msg )
    sys.exit()


def main(archive, archivename, logfile, logid, workdir, merge, rm_header=0, concat=''):

    # create a temporary repository
    tmp_dir = tempfile.mkdtemp(dir=os.environ['GALAXY_HOME']+'/database/tmp/')

    #open log file
    mylog = open(logfile, "w");

    is_gzipped, is_gzvalid = check_gzip( archive )
    is_bzipped, is_bzvalid = check_bz2( archive )

    # extract all files in a temp directory
    # test if is a zip file
    if check_zip( archive ):
	with zipfile.ZipFile(archive, 'r') as myarchive:
		myarchive.extractall(tmp_dir)

    # test if is a tar file
    elif tarfile.is_tarfile( archive ):
	mytarfile=tarfile.TarFile.open(archive)
	mytarfile.extractall(tmp_dir)
	mytarfile.close()
 
    # test if is a gzip file
    elif is_gzipped and is_gzvalid :
	mygzfile = gzip.open(archive, 'rb')

	myungzippedfile = open (tmp_dir+"/"+os.path.splitext(os.path.basename(archivename))[0], 'wb', 2**20)
	for i in iter(lambda: mygzfile.read(2**20), ''):
		myungzippedfile.write(i) 

	myungzippedfile.close()
	mygzfile.close()

    elif is_bzipped and is_bzvalid:
        mybzfile = bz2.BZ2File(archive, 'rb')

        myunbzippedfile = open (tmp_dir+"/"+os.path.splitext(os.path.basename(archivename))[0], 'wb', 2**20)
        for i in iter(lambda: mybzfile.read(2**20), ''):
                myunbzippedfile.write(i)

        myunbzippedfile.close()
        mybzfile.close()

		
    # test if merge is enable
    if merge == "true":
	mylog.write("Merge option is enabled with "+str(rm_header)+" lines to deleted\n\n")
	myfinalfile = open(concat, "w");
	for myfile in listdirectory(tmp_dir):
		myopenfile = open(myfile, "r")
		nblinesremove=0
		mylog.write(myfile+" is extracted from the archive and is added into the result file\n")
		for line in myopenfile:	

			#if not equal, don't write	
			if int(rm_header) != nblinesremove:
				nblinesremove=nblinesremove+1
			else:
				# write the line into the final file			
				myfinalfile.write(line)
	
	myfinalfile.close()

    else:
	# if merge is disable
        mylog.write("Merge option is disabled\n\n")
    	my_files={}

   	# move all files (recursively) in the working dir
    	for myfile in listdirectory(tmp_dir):
		# /opt/path/sub dir/my_file.txt -> sub\ dir/my.file.txt
		myfile = myfile.split(tmp_dir)[1].split('/',1)[1].replace(" ", "\ ")
		myfileclean = os.path.basename(myfile).replace("_", ".")

		# verify and add file name into a tab
		while myfileclean in my_files:
			myfileclean=os.path.splitext(myfileclean)[0]+'1'+os.path.splitext(myfileclean)[1]

		my_files[myfileclean]='ok'

		# if rename file
		if ( myfileclean != myfile):
			mylog.write(myfile+" is extracted from the archive and is renamed : "+myfileclean+"\n")
		else:
    			mylog.write(myfile+" is extracted from the archive\n")

		fileext = os.path.splitext(myfile)[1].replace(".", "")

		# if no extension
		if fileext == '':
			fileext="txt"

		if fileext == 'fa':
			fileext="fasta"

		if fileext == 'fq':
			fileext="fastq"

		# move files in working dir
		os.system("mv "+tmp_dir+"/"+myfile+" "+workdir+"/primary_"+logid+"_"+myfileclean+"_visible_"+fileext)

   	mylog.write("\n\nPlease refresh your history if all files are not present\n")
    	mylog.close()
	

    #clean up temp files
    shutil.rmtree( tmp_dir, ignore_errors=True )


# parse the directory and return files path (in a tab)
def listdirectory(path):
	myfile=[]
	l = glob.glob(path+'/*') 
    	for i in l:
		# if directory
        	if os.path.isdir(i): 
			myfile.extend(listdirectory(i))
		# else put the file in the tab      	  	
		else:
			myfile.append(i)
    	return myfile


if __name__=="__main__": main(*sys.argv[1:])
