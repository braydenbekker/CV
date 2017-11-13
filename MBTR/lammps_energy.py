'''
Authors Chandramouli Nyshadham, Brayden Bekker

Generate the ground state energies of fcc, bcc. or hcp based structures 
from a list of binary alloys (i.e. Cu Ta) CuTa.

The structures are generated using enumlib. 

We will be using Quippy and lammpslib package. 

The interatomic potentials for the binary systems were downloaded from NIST interatomic 
potentials database.

Use of docker, quippy and lammpslib package was taught to me by Conrad Rosenbrock.
'''
import quippy
from lammpslib import LAMMPSlib
#from aflow import *
import os
from os.path import isfile, join
import numpy as np
import collections
import matplotlib.pyplot as plt

# loading the interatomic potential.

header = ["units metal",
          "dimension 3",
          "boundary p p p",
          "atom_style atomic",
          "atom_modify map array"]
cmds = ["pair_style adp",
        "pair_coeff * * /fslhome/bbekker/MBTR/eam_potentials/Al-Ti-eam.alloy Al Ti", # check this
        "mass 1 26.981539", #Al
        "mass 2 47.867", #Ti
        "neighbor 2.0 bin",
        "neigh_modify delay 10 check yes"]
pot = LAMMPSlib(lmpcmds=cmds, atom_types={"Al": 1, "Ti": 2}, lammps_header=header, keep_alive=True)

pathToDirectories=os.getcwd()
folders=[pathToDirectories+str(i) for i in ['/Structures/fcc/','/Structures/bcc/','/Structures/hcp/']]
folders

def readAllVASPFiles(paths,potential,numStructToRead=0,AtomicNum=[13,22]):
    """
    This fuction reads in all the crystal structures in the path to folder given. All the input files should be in
    vasp format. 
    
    The function returns atom numbers, positions(cartesian coordinates), and total energies computed using the
    potential.
    
    Arguments:
        Input: 
            path (string): path to the folder where all strucutres (vasp files) are present
            numStructToRead (int): if FALSE all files in the folder are read,else only 'numStructToRead' are read.
            potential: potential for computing the total energies.
            AtomicNum: Array containing elements atomic numbers.
            
        output:
            Returns 
            
            z (list): atomic nuclear charge
            pos (list): positions  of all atoms in a crystal structure. (cartesian) 
            lattice (list): lattice vectors of the unit cell. 
            totEnePerAtom (list): total energy per atom of crystal structure computed using potential.
            structInfo (list): contains the structure information (file number).
            conc : concentration of element A in the binary AB (listed in alphabetical order) 
            
    Note: All files in the folder should be named in the format, 'vasp.1', 'vasp.3', ... etc. which is what enumlib 
            generates.
    """
    
    AtomicNumA=AtomicNum[0]
    AtomicNumB=AtomicNum[1]

    res_z=[]
    res_pos=[]
    res_lattice=[]
    res_totEnePerAtom=[]
    res_structInfo=[]
    res_conc=[]
    
    res={}
  
    for path in paths:
        if numStructToRead:
            totalFiles=numStructToRead+1

        else:
            totalFiles=len(os.listdir(path))

        for i in range(1,totalFiles):

            inpFile=join(path,'vasp.'+str(i)) # get the file name
            a = quippy.Atoms(inpFile, format="vasp") # read the input file as quippy object 
            a.set_calculator(potential) # set the potential 

            res_z.append(list(a.z)) # get nuclear charge
            res_pos.append(a.positions.tolist()) # get positions in cartesian
            res_lattice.append(a.cell.tolist())  # get lattice vectors

            totalNumOfAtoms=len(a.z) # total number of atoms
            #print a.get_total_energy()
            #res_totEnePerAtom.append(a.get_total_energy()/float(len(a.z))) # compute total energy per atom

            conc=a.z.tolist().count(AtomicNumB)/float(len(a.z)) # get the concentration.
            res_conc.append(conc)

            res_structInfo.append(path[-4::]+str(i)) # structure information a.k.a file name.

    res={'z':res_z,'pos':res_pos,'lattice':res_lattice,'totEnePerAtom':res_totEnePerAtom,'structInfo':res_structInfo,'conc':res_conc}

    return res

def computeEnthalpy(data):
    """
    
    Input:
    data : Dictionary generated using 
    readAllVASPFiles(paths,potential,numStructToRead=0,AtomicNum=[29,73]) function.'
    
    Output:
    enthalpy added to dictionary 'data'.
    
    """
    
    res_enthalpy=[]
    
    # compute the lowest total energies of pure elements.
    
    # element A with concentration 0.0
    index=[i for i,x in enumerate(data['conc']) if x == 0.0]
    fileNum=data['totEnePerAtom'].index(np.array(data['totEnePerAtom'])[index].min())
    ene_A=data['totEnePerAtom'][fileNum]
    print "element A", data['structInfo'][fileNum], data['totEnePerAtom'][fileNum]
    
     # element B with concentration 1.0
    index=[i for i,x in enumerate(data['conc']) if x == 1.0]
    fileNum=data['totEnePerAtom'].index(np.array(data['totEnePerAtom'])[index].min())
    ene_B=data['totEnePerAtom'][fileNum] 
    print "element B", data['structInfo'][fileNum], data['totEnePerAtom'][fileNum]
    
    # compute formation enthalpy per atom.
    # formulae: \Delta Hf = total_energy_per_atom_of_structure - sum_of_total_energies_of pure_elements
    
    for i in range(len(data['conc'])):
        enthalpy = data['totEnePerAtom'][i] - (1-data['conc'][i])*float(ene_A)  - (data['conc'][i])*ene_B
        res_enthalpy.append(enthalpy)
    
    return res_enthalpy

dat=readAllVASPFiles(folders,potential=pot)

dat.update({'enthalpy':computeEnthalpy(dat)})

# Seperate FCC, BCC, HCP

def plotConvexHull(data):
    """
    Plots the convex hull.
    """
    
    fccNum=len(filter(lambda x:'FCC' in x, data['structInfo']))
    bccNum=len(filter(lambda x:'BCC' in x, data['structInfo']))
    hcpNum=len(filter(lambda x:'HCP' in x, data['structInfo']))
            
    # plotting the convex hull
    plt.figure()
    plt.scatter(data['conc'][0:fccNum-1],data['enthalpy'][0:fccNum-1],label='FCC')
    plt.scatter(data['conc'][fccNum:fccNum+bccNum-1],data['enthalpy'][fccNum:fccNum+bccNum-1],label='BCC')
    plt.scatter(data['conc'][-hcpNum:],data['enthalpy'][-hcpNum:],label='HCP')
    plt.xlim(0,1)
    plt.legend()
    plt.show()

plotConvexHull(dat)
