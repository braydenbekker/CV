'''
Author Brayden Bekker

getStructures.py uses the enumlib code created by the MSG group at BYU, https://github.com/msg-byu/enumlib. enumlib is 
designed to compute the derivative superstructures for a system of elements and a given lattice. 
Input:
    - path to makeStr.py and enum.x from enumlib.
    - precomputed struct_enum.out files for FCC, BCC, and HCP through 7140 structures.
Output:
    - a directory (e.g. Structures_AlTi) whith sub directories, fcc, bcc, hcp
    - each subdirectory contains 2500 vasp POSCAR files with atom concentrations and positions.
    (the first 500 derivative superstructures and 2000 additional structures from higher atom cells 10-12)
    
File Structure:
cd path/to/program/files/
ls
Structures/ struct_enum.in.fcc, struct_enum.in.bcc, struct_enum.in.hcp
enum_(fcc,bcc,hcp)/struct_enum.out
makeStr.py
./enum.x
'''
import shutil
import numpy as np
import os
import random
import sys
def GenEnumOut(structure):
    '''
    Generate the struct_enum.out files if this has not been done previously.
    '''
    input_file=('struc_enum.in.%s'%(structure)) #FCC, BCC, HCP input structure                                                              
    directory = "/fslhome/bbekker/MBTR" #file path to enumlib files                                                                         
    os.chdir(directory)
    #run the enum.x fortran script wiht struct_enum.in.sys as input to generate the struct_enum.out file.                                    
    enum="./enum.x %s" %(input_file)
    os.system(enum)
    return

def GenVaspFiles(system, structure, num, Species): #(Al Ti, fcc, 500, Al_Ti)
    '''
    Generate 2500 vasp POSCAR files for a given system of elements and structure
    '''
    directory = ("/fslhome/bbekker/MBTR/Structures/%s"%(structure)) #file path to enumlib files
    os.chdir(directory)

    os.system(("cp ~/MBTR/enum_%s/struct_enum.out %s") %(structure, directory))
    
    makestructures="python ~/MBTR/makeStr.py 1211 7140 2000 -species %s" %(system)
    os.system(makestructures)

    makestructures="python ~/MBTR/makeStr.py 1 %s -species %s" %(num, system)
    os.system(makestructures)
    
    os.system("rm struct_enum.out")
    
def main():
    '''
    setup directories for computed vasp POSCAR files based on command line args.
    '''
    num = 500
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
