import Bio
from Bio.PDB import *
import numpy as np
from itertools import repeat

# Input the appropriate filename (including .pdb extension) for the structure you want to look at

parser = PDBParser(PERMISSIVE=1)

def charge_site_calc(file, charge_state):
    
    filename = file
    charge = charge_state
    structure_id = filename.replace('.pdb','')
    structure = parser.get_structure(structure_id, filename)

    model = structure[0]
    resname_list = ['LYS', 'ARG', 'ASP', 'GLU']
    resposi_list = []
    prot_list = []
    deprot_list = []
    charge_list = []
    
    lys_counter = 0
    arg_counter = 0
    asp_counter = 0
    glu_counter = 0
    his_counter = 0
    h2o_counter = 0
    het_counter = 0
    h2o_counter = 0
    het_counter = 0

    for model in structure:
        for chain in model:     
            for residue in chain:
                if residue.resname == 'LYS':
                    lys_counter = lys_counter + 1
                if residue.resname == 'ARG':
                    arg_counter = arg_counter + 1
                if residue.resname == 'ASP':
                    asp_counter = asp_counter + 1
                if residue.resname == 'GLU':
                    glu_counter = glu_counter + 1
                if residue.resname == 'HIS':
                    his_counter = his_counter + 1

            lys_string = ('LYS: ' + str(lys_counter))
            arg_string = ('ARG: ' + str(arg_counter))
            asp_string = ('ASP: ' + str(asp_counter))
            glu_string = ('GLU: ' + str(glu_counter))
            his_string = ('HIS: ' + str(his_counter))

            total_counter = (lys_counter + arg_counter + asp_counter + glu_counter + his_counter)

            print('Number of residues in chain ' + str(chain.get_id()) +' available to be protonated:')
            print(lys_string)
            print(arg_string)
            print(asp_string)
            print(glu_string)
            print(his_string)
            print('Total: ' + str(total_counter))

            # Creating lists for the pdb2gmx_list.txt for proton assignment
            # First, create list of ionisable sites for each chain
            resposi_list.extend(repeat('LYS', lys_counter))
            resposi_list.extend(repeat('ARG', arg_counter))
            resposi_list.extend(repeat('ASP', asp_counter))
            resposi_list.extend(repeat('GLU', glu_counter))
            resposi_list.extend(repeat('HIS', his_counter))
            resposi_list.append('NTM')
            resposi_list.append('CTM')
            
            # Create list turning ionisation sites into protonation values
            
            for item in resposi_list:
                if item in resname_list:
                    prot_list.append(1)
                    deprot_list.append(0)
                if item == 'HIS':
                    prot_list.append(1)
                    deprot_list.append(2)
                if item == 'NTM':
                    prot_list.append(0)
                    deprot_list.append(2)
                if item == 'CTM':
                    prot_list.append(2)
                    deprot_list.append(0)
            
            # Reset residue counters for new chain
            
            lys_counter = 0
            arg_counter = 0
            asp_counter = 0
            glu_counter = 0
            his_counter = 0
            h2o_counter = 0
            het_counter = 0
            
    # Checking heteroatoms and waters
    
    for residue in chain.get_list():
        residue_id = residue.get_id()
        hetfield = residue_id[0]
        if hetfield[0] == 'H':
            het_counter = het_counter + 1
        if hetfield[0] == 'W':
            h2o_counter = h2o_counter + 1

    print('Number of heteroatoms: ' + str(het_counter))
    print('Number of waters: ' + str(h2o_counter))

    if het_counter or h2o_counter != 0:
        print('Heteroatoms and water molecules will need to be removed prior to running the algorithm')
            
    # Creating list of protonated sites to get specific charge state for protein,
    # residues chosen in order, randomised positions not required as the charges
    # will spread around the surface once the algorithm starts.
    
    for item in resposi_list:
        if charge >= 1:
            if item in resname_list:
                charge_list.append(1)
                charge = (charge - 1)
            if item == 'HIS':
                charge_list.append(1)
                charge = (charge - 1)
            if item == 'NTM':
                charge_list.append(1)
                charge = (charge - 1)
            if item == 'CTM':
                charge_list.append(1)
                charge = (charge - 1)         
        elif charge == 0:
            if item in resname_list:
                charge_list.append(0)
            if item == 'HIS':
                charge_list.append(2)
            if item == 'NTM':
                charge_list.append(2)
            if item == 'CTM':
                charge_list.append(0)
    
    # Write out protonation and deprotonation lists for use with the algorithm

    with open(structure_id + '_prot_list.txt', 'w') as file:
        for item in prot_list:
            file.write(str(item) +'\n')
    print('Protonated list written')

    with open(structure_id + '_deprot_list.txt', 'w') as file:
        for item in deprot_list:
            file.write(str(item) +'\n')
    print('Deprotonated list written')
            
    with open(structure_id + '_' + str(charge_state) + '_list.txt', 'w') as file:
        for item in charge_list:
            file.write(str(item) +'\n')
    print('Protonation list for a +' + str(charge_state) + ' protein written')
            
    return