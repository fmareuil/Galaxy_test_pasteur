#!/usr/bin/env python
# library creation, links directory creation, outputs directory creation and chown output data
import os, sys, string, urllib2, smtplib, stat
sys.path.insert( 0, os.path.dirname( __file__ ) )
from common import display, get, post, submit
from pwd import getpwnam, getpwuid
from email.mime.text import MIMEText
from collections import namedtuple

GALAXY_EMAIL = "Galaxy_pasteur@pasteur.fr"
#ADMINS_EMAILS = ["fmareuil@pasteur.fr"]
ADMINS_EMAILS = ["fmareuil@pasteur.fr","odoppelt@pasteur.fr","screno@pasteur.fr"]
COMMASPACE = ", "


def main():
    try:
        key = sys.argv[1]
        base_url = sys.argv[2]
        base_link = sys.argv[3]
        base_output = sys.argv[4]
	try:
            gid = int(sys.argv[5])
        except ValueError:
            print >> sys.stderr, "Error, gid must be an integer"
            sys.exit( 1 )
    except IndexError:
        print >> sys.stderr, 'usage: %s key base_url base_link base_output gid' % os.path.basename( sys.argv[0] )
        sys.exit( 1 )  
    return create_libraries_folders( key, base_url, base_link, base_output, gid)

def sendmail(msg,subject):
    # Create a text/plain message
    msg = MIMEText(msg)
    # me == the sender's email address
    # you == the recipient's email address
    msg['Subject'] = subject
    msg['From'] = GALAXY_EMAIL
    msg['To'] = COMMASPACE.join(ADMINS_EMAILS)

    # Send the message via our own SMTP server, but don't include the
    # envelope header.
    s = smtplib.SMTP('smtp.pasteur.fr')
    s.sendmail(GALAXY_EMAIL, ADMINS_EMAILS, msg.as_string())
    s.quit()   


def create_libraries_folders(key, base_url, base_link, base_output, gid):
    try:
        folders_output = os.listdir(base_output)
    except OSError, e:
        print >> sys.stderr, ("Error, %s doesn't exist" % base_output)
        sys.exit( 1 )

    libcreation='Galaxy library creation'
    
    for j in folders_output:
                
        directory = os.path.join(base_output, j)
        user_uid = os.stat(directory).st_uid
        # exception is caught and an email is sent, if a uid is no longer in LDAP
        try:

           
            user_gid = getpwuid(user_uid).pw_gid
            files = os.listdir(directory)
           
        #    print ("Debug cron pb, user: %s - user_uid: %s " % (j,user_uid))
            for r in files:
                try:
                    file = os.path.join(directory, r)
                    stat_file = os.stat(file)
                    uid_file = stat_file.st_uid
                    if user_uid != uid_file:
                    #pass
                # TEMPORARY FIX exportdata 
                # Rights are modified during the file export therefor it fails
                    #st=os.stat(file)
                    
                    #os.chmod(file, st.st_mode | stat.S_IROTH)
            
                        os.chown(file,user_uid,user_gid)
                except Exception, e:
                    sendmail("Error during chown of the file in the output directory of %s user : %s" % (j, e),libcreation)
                    sys.exit( 1 )
        except Exception, e:
            sendmail("LDAP problem for user: %s - user_uid: %s " % (j,user_uid),'Galaxy: LDAP user problem')
            #sendmail ("LDAP problem for user: %s - user_uid: %s " % (ldap_pbs_users,ldap_pbs_uid), 'Galaxy: LDAP user problem')
            
#sys.exit(1)
    try:
        users = display(key, base_url+"users", return_formatted=False)
    except Exception, e:
        print >> sys.stderr, ("Error during the display API galaxy command : %s" % e)
        sys.exit( 1 )
    try:
        folders = os.listdir(base_link)
    except OSError, e:
        print >> sys.stderr, ("Error, %s doesn't exist" % base_link)
        sys.exit( 1 )
    for i in users:
        if not i["email"] in folders:
            data = {}
            data[ 'name' ] = i["email"]
            try:
                folder = os.path.join(base_link,i["email"])
                os.mkdir(folder, 02750)
                login = string.split(i["email"],"@")[0]
                uid = getpwnam(login).pw_uid
                os.chown(folder,uid,gid)
                data[ 'description' ] = "Library of %s" % login
                url_libraries = os.path.join(base_url,"libraries")
                lib = submit( key, url_libraries, data, return_formatted=False)
		datapermission = {}
                datapermission['LIBRARY_ACCESS_in'] = i['id']
                datapermission['LIBRARY_MODIFY_in'] = i['id']
                datapermission['LIBRARY_MANAGE_in'] = i['id']
                datapermission['LIBRARY_ADD_in'] = i['id']
                url_permission = os.path.join(base_url,"libraries",lib['id'],"permissions")
                submit( key, url_permission, datapermission )
                sendmail("Creation of the library of %s user" % (i["email"]),libcreation)
            except Exception, e:
                sendmail(("Error during the creation of the directory and library of %s user : %s" % (i["email"], e)),libcreation)
                sys.exit( 1 )
        if not i["email"] in folders_output:
            try:
                folder_output = os.path.join(base_output,i["email"])
                os.umask(0)
                os.mkdir(folder_output, 02770)
                login = string.split(i["email"],"@")[0]
                uid = getpwnam(login).pw_uid
                os.chown(folder_output,uid,gid)
            except Exception, e:
                sendmail("Error during the creation of the output directory of %s user : %s" % (i["email"], e),libcreation)
                sys.exit( 1 )

if __name__ == "__main__":
    main()


