################################
#                                                     #
#  export.py    version 2.2                 #
#                                                     #
#                                                     #
#                                                     #
#                                                     #
#  Galaxy Team copyright                 #
#  CIB, Institut Pasteur                     #
###############################
#!/usr/bin/env python


import optparse, re, sys, os, shutil

class directory:
  '''
  Class that automatically traverses directories
  and builds a tree with size info
  '''
  def __init__(self, path, parent=None):

    if path[-1] != '/':
        # Add trailing /
        self.path = path + '/'
    else:
        self.path = path
    self.size = 4096
    self.parent = parent
    self.children = []
    self.errors = []
    for i in os.listdir(self.path):
        try:
            self.size += os.lstat(self.path + i).st_size
            if os.path.isdir(self.path + i) and not os.path.islink(self.path + i):
                a = directory(self.path + i, self)
                self.size += a.size
                self.children.append(a)
        except OSError:
            self.errors.append(path + i)
    self.size = self.size/1024
	    
	    
logout = open("output_export.log" , "w")

def foo_callback(option, opt, value, parser):
    setattr(parser.values, option.dest, value.split())

def name_test(obj_opt_name):
    regex = re.compile("[A-Za-z0-9_\-.]")
    for i in obj_opt_name:
        if regex.sub('', i) != "":
	    print >> sys.stderr, "wrong character in one of the filename"
	    sys.exit( 1 )
            return False
	else:
	    return True

def copy_files(obj_opt):
    for i in range(len(obj_opt.name)):
        if os.path.exists(obj_opt.export_dir+"/"+obj_opt.name[i]):
            print >> sys.stderr, "%s/%s already exists: delete file or change the name" % (obj_opt.export_dir,obj_opt.name[i])
	    sys.exit( 1 )
        st = directory(obj_opt.export_dir)
        used = st.size    
        if used >= 104857600:
            print >> sys.stderr, "%s is too large: move or delete files \nif you do not have permission, please contact an admin" % (obj_opt.export_dir)
            sys.exit( 1 )
        else:
            logout.write("command : cp %s %s/%s\n" % (obj_opt.input_name[i], obj_opt.export_dir, obj_opt.name[i]))
            #shutil.copyfile(obj_opt.input[i], obj_opt.export_dir+"/"+obj_opt.name[i])
	    outputtmp = obj_opt.export_dir_tmp+"/"+obj_opt.user+"_"+obj_opt.name[i]
            shutil.copyfile(obj_opt.input[i], outputtmp)
            shutil.move(outputtmp, obj_opt.export_dir+"/"+obj_opt.name[i])
            os.chmod(obj_opt.export_dir+"/"+obj_opt.name[i],0660)
            logout.write("%s \nhave been copied in \n%s \nand named \n%s \nWARNING: this tool duplicates the data, remember to delete redundant data.\n" % (obj_opt.input_name[i],obj_opt.export_dir,obj_opt.name[i]))	     
    logout.close()

if __name__ == '__main__':
    op = optparse.OptionParser()
    op.add_option('-i', '--input', type='string', action='callback',callback=foo_callback)
    op.add_option('-f', '--input_name', type='string', action='callback',callback=foo_callback)
    op.add_option('-n', '--name', type='string', action='callback',callback=foo_callback)
    op.add_option('-u', '--user', default=None)
    op.add_option('-e', '--export_dir', default=None)
    op.add_option('-o', '--output', default=None)
    #print help(op)
    opts, args = op.parse_args()
    opts.export_dir_tmp = opts.export_dir+"/../outputs/.tmp"
    opts.export_dir = opts.export_dir+"/../outputs/"+opts.user
    #print opts.name  
    if name_test(opts.name):
        try:
            copy_files(opts)
        except Exception, e:
            print >> sys.stderr, e
            sys.exit( 1 )
   


 
