#!/usr/bin/python
import string
import os
import sys
import argparse

class Base():
	
    def verif(self):
        if os.path.exists(self.path_loc):
            fileloc = open(self.path_loc, 'r')
	    for i in fileloc:
	        line = string.split(i)
		if len(line) > 1:
	            if line[self.id] == self.name:
	                fileloc.close()
	                return 1    
    	    fileloc.close()
	    return 0
	else:
	    return 0

    def upgrade(self):
        path_loctmp = self.path_loc + "tmp"
	fileloctmp = open(path_loctmp ,"w")
        fileloc = open(self.path_loc, 'r')
	for i in fileloc:
	    line = string.split(i)
	    if len(line) > 0:
	        if self.name == line[self.id]:
	            fileloctmp.write(self.structure)
	        else:
	            fileloctmp.write(i)
	fileloctmp.close()
	fileloc.close()
	os.rename(path_loctmp, self.path_loc)
		
    def add(self):
        alreadyexist = self.verif()
	if alreadyexist:
	    self.upgrade()
	    print  >> sys.stdout, "Path %s for %s updated" % (self.path_genome, self.name)
	else:
            file = open(self.path_loc, 'a')
	    file.write(self.structure)
            file.close()
	    print  >> sys.stdout, "Path %s for %s added" % (self.path_genome, self.name)

    def remove(self):
        pass

	
class GatkPicardBwaBowtie(Base):
    def __init__(self, name, path_loc, path_genome):
        self.path_loc = path_loc
	self.path_genome = path_genome
	self.name = name
	self.structure = "%s\t%s\t%s\t%s\n" % (self.name, self.name, self.name, self.path_genome)
	self.id = 0

class Soap(Base):
    def __init__(self, name, path_loc, path_genome):
        self.path_loc = path_loc
	self.path_genome = path_genome
	self.name = name
	self.structure = "%s\t%s\n" % (self.path_genome, self.name)
	self.id = 1


##########################################
##################MAIN###################
##########################################
	    
def path(string):
    if os.path.exists(string):
        return string
    else:
        print  >> sys.stderr, "Error, %s does not exist" % (string)
	sys.exit( 1 )
 
def file_conf(args):
     conf_file = open(args.f, "r")
     for j in conf_file:
         args.n, args.x, args.g, args.L, args.l = None, None, None, None, None
         line = string.split(j)
	 if len(line) > 3:
	     args.n = line[0]
             args.x = line[1]
	     args.g = line[2]
	     path(line[2])
	     if len(line) == 5:
	         if line[4] == "create":
	             args.L = line[3]
	     else:
	         args.l = line[3]
	         path(line[3])
	 else:
	     print  >> sys.stderr, "Error, %s is not a valid configuration line" % (j)
	     sys.exit( 1 )
	 manager(args)
     conf_file.close()
        
def manager(arguments): 
    if arguments.L:
        arguments.l = arguments.L
	
	
    if arguments.x == "bwa" or arguments.x == "bowtie" or arguments.x == "picard"  or arguments.x == "gatk" or arguments.x == "gatk2" or arguments.x == "bowtie2":
        objectindex = GatkPicardBwaBowtie(arguments.n, arguments.l, arguments.g)
        objectindex.add()
    
    elif  arguments.x == "soap":
        objectindex = Soap(arguments.n, arguments.l, arguments.g)
	objectindex.add()
	
    else:
        print  >> sys.stderr, "Error, %s is not a valid indexe" % (arguments.x)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='file configuration', type=str)
    parser.add_argument('-x', help='type of indexes, bwa, picard, bowtie, bowtie2, soap', type=str)
    parser.add_argument('-n', help='name of genome', type=str)
    parser.add_argument('-l', help='path of the loc file', type=path)
    parser.add_argument('-L', help='path of the new loc file', type=str)
    parser.add_argument('-g', help='path of the genome indexes', type=path)
    args = parser.parse_args()
    if args.f:
        file_conf(args)
    else:
        manager(args)
