#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 19:37:05 2018

@author: Haotian.Teng
"""
def make_markdown_table(array):

    """Generate a markdown format string given a python nested list. 
    Args: Python list with rows of table as lists
               First element as header. 
        Output: String to put into a .md file 
        
    Ex Input: 
        [["Name", "Age", "Height"],
         ["Jake", 20, 5'10],
         ["Mary", 21, 5'7]] 
    """


    markdown = "  \n" + str("|")

    for e in array[0]:
        to_add = str(e) + str("|")
        markdown += to_add
    markdown += "  \n"

    markdown += '|'
    for i in range(len(array[0])):
        markdown += str("-|")
    markdown += "  \n"

    for entry in array[1:]:
        markdown += str("|")
        for e in entry:
            to_add = str(e) + str("|")
            markdown += to_add
        markdown += "  \n"

    return markdown + "  \n"

def mark_down_print_protein(protein):
    """
    This functions is used to generate markdown format report string of a given protein 
    Args:
        protein: A Protien class.
    Return:
        md_string:String content the markdown content.
    """
    content = ''
    content += '# '+ protein.pdb_id + '  \n'
    for idx,chain in enumerate(protein.chains):
        content += '## chain'+ chr(idx+65) + '  \n'
        content += '### Table 1'+chr(idx+65) + ' Bond Length and Supplemental Bond Angle Theta  \n'
        coor_table = [['Index','Atom', 'Bond_Length','Theta','x','y','z']]
        for atom_idx,coor in enumerate(chain.coor):
            bond = []
            bond.append(atom_idx)
            bond.append(chain.atoms[atom_idx])
            if atom_idx == 0:
                bond.append('nan')
            else:
                bond.append(chain.bond[atom_idx-1])
            if (atom_idx == 0) or (atom_idx == (len(chain.atoms)-1)):
                bond.append('nan')
            else:
                bond.append(chain.theta[atom_idx-1])
            bond+=chain.coor[atom_idx]
            coor_table.append(bond)
        content += make_markdown_table(coor_table)
        content += '### Table 2'+chr(idx+65) + ' The Three Dihedrals  \n'
        dihedrals_table = [['Index','Residue','Phi','Psi','Omega']]
        for aa_idx, aa in enumerate(chain.aas):
            d = []
            d+= [aa_idx,aa]
            d+= list(chain.dihedrals[aa_idx])
            dihedrals_table.append(d)
        content += make_markdown_table(dihedrals_table)
    return content