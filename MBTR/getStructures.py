'''
Author Brayden Bekker

Get the Structure and Save to a file
'''
import shutil
import numpy as np
import os
import random
import sys
def GenEnumOut(structure):
    input_file=('struc_enum.in.%s'%(structure)) #FCC, BCC, HCP input structure                                                              
    directory = "/fslhome/bbekker/MBTR" #file path to enumlib files                                                                         
    os.chdir(directory)
    #run the enum.x fortran script wiht struct_enum.in.sys as input to generate the struct_enum.out file.                                    
    enum="./enum.x %s" %(input_file)
    os.system(enum)
    return

def GenVaspFiles(system, structure, num, Species):
    directory = ("/fslhome/bbekker/MBTR/Structures/%s"%(structure)) #file path to enumlib files
    os.chdir(directory)

    os.system(("cp ~/MBTR/enum_%s/struct_enum.out %s") %(structure, directory))
    
    #makestructures="python ~/MBTR/makeStr.py 1211 7140 2000 -species %s" %(system)
    #os.system(makestructures)

    #makestructures="python ~/MBTR/makeStr.py 1 %s -species %s" %(num, system)
    #os.system(makestructures)

    makestructures="python ~/MBTR/makeStr.py 1211 7140 5 -species %s" %(system)                                       
    os.system(makestructures)

    os.system("rm struct_enum.out")
    
def main():
    num = 500 #change later when working to 20000
    structure = np.array(["fcc","bcc","hcp"])
    
    first=sys.argv[1]
    second=sys.argv[2]
    species="%s_%s"%(first, second)
    system="%s %s"%(first, second)

    os.mkdir("Structures")
    os.mkdir("./Structures/fcc")
    os.mkdir("./Structures/bcc")
    os.mkdir("./Structures/hcp")

    for i in range(3):
        GenVaspFiles(system, structure[i], num, species)

    # os.chdir("~/MBTR")
    #filename=("cp -rf Structures/ Structures_%s%s" %(first,second))
    #os.system(filename)
    #os.system("rm Structures/")
    
main()