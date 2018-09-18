#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 21:30:25 2018

@author: Haotian.Teng
"""
import utils.geometry_op as g_op
from utils.md_print import mark_down_print_protein
from collections import defaultdict
import numpy as np
import os
from matplotlib import pyplot as plt
COLOR_PLATTER = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
class Chain(object):
    def __init__(self):
        self.aas = list()
        self.atoms = list()
        self.aa_index = list()
        self.coor = list()
        self.bond = list()
        self.theta = list()
        self.dihedrals = list()
        """
        Attrs:
            aas: list of amino acids sequences, aa is represented as three capital letter.
            atoms: the list of atoms in the chain, from N to C.
            aa_index: the index of residue each atom belongs to.
            coor: Nx3 list, contain the x,y,z-coordinates of each atom.
            bond: N-1 list, contain the lenght of bond between the adjacent atoms.
            theta: N-2 list, contain the angle of the current atom and its adjacent atoms.
            dihedrals: Rx3 list, R is the number of the residue, given the dihedrals
                of each residue in the [phi,psi,omega] order.
        """
    def add_aas(self,aa_segments):
        self.aas+=aa_segments
    def add_atom(self,atom_info):
        """
        Args:
            atom_info: a list of str, the ATOM split line in pdb file.
        """
        if atom_info[0] != 'ATOM':
            raise ValueError("Input must be the ATOM line of pdb file.")
        atom = atom_info[2][0:4]
        self.atoms.append(atom)
        if len(atom_info[3]) == 1:
            self.aa_index.append(int(atom_info[4])-1)
            self.coor.append([float(x) for x in atom_info[5:8]])
        elif len(atom_info[4])==1:
            self.aa_index.append(int(atom_info[5])-1)
            self.coor.append([float(x) for x in atom_info[6:9]])
    def _map_Cartesian(self):
        r_dd = lambda: defaultdict(r_dd)
        self.aa_index_dict = r_dd()
        for idx, aa_idx in enumerate(self.aa_index):
            self.aa_index_dict[aa_idx][self.atoms[idx]]=idx
    def parse(self):
        """
        This will calculate the bond length, theta and Dihedrals based on current state of the Chain
        """
        self.bond = list()
        self.theta = list()
        self.Dihedrals = list()
        self._map_Cartesian()
        for idx,c in enumerate(self.coor[1:-1]):
            self.bond.append(g_op.get_len(self.coor[idx],c))
            self.theta.append(g_op.get_angle(self.coor[idx],c,self.coor[idx+2]))
        for idx,_ in enumerate(self.aas):
            N = self.coor[self.aa_index_dict[idx]['N']]
            CA = self.coor[self.aa_index_dict[idx]['CA']]
            C = self.coor[self.aa_index_dict[idx]['C']]
            if idx == 0:
                phi = np.nan
            else:
                C0 = self.coor[self.aa_index_dict[idx-1]['C']]
                phi = g_op.get_diherdral([C0,N,CA,C])
            if idx == len(self.aas)-1:
                psi = np.nan
                omega = np.nan
            else:
                N1 = self.coor[self.aa_index_dict[idx+1]['N']]
                CA1 = self.coor[self.aa_index_dict[idx+1]['CA']]
                psi = g_op.get_diherdral([N,CA,C,N1])
                omega = g_op.get_diherdral([CA,C,N1,CA1])
            self.dihedrals.append(np.asarray([phi,psi,omega]))
                
        self.bond.append(g_op.get_len(c,self.coor[idx+2]))
class Protein(object):
    def __init__(self):
        self.pdb_id = None
        self.chains = []
    def add_chain(self,chain_instance):
        if not isinstance(chain_instance,Chain):
            raise TypeError("Expect a chain instance but a %s is given"%(type(chain_instance)))
        self.chains.append(chain_instance)
    def get_chain(self,idx):
        if type(idx) is int:
            return self.chains[idx]
        elif type(idx) is str:
            if len(idx) > 1:
                raise IndexError("Chain index must be either integer or character, but string %s found."%(idx))
            else:
                idx = ord(idx)
            if (idx>122) or (idx<65):
                raise IndexError("Chain index must be character from A-Z or a-z, but character %s found."%(idx))
            if idx>96:
                return self.chains[idx - ord('a')]
            else:
                return self.chains[idx - ord('A')]
    def read_pdb(self,pdb_file):
        with open(pdb_file,'r') as f:
            for line in f:
                split_line = line.split()
                if split_line[0] == 'HEADER':
                    self.pdb_id = split_line[-1]
                if split_line[0] == 'SEQRES':
                    if int(split_line[1]) == 1:
                        temp_chain = Chain()
                        temp_chain.add_aas(split_line[4:])
                        self.add_chain(temp_chain)
                    else:
                        self.get_chain(split_line[2]).add_aas(split_line[4:])
                if line.startswith('ATOM'):
                    try:
                        chain = self.get_chain(split_line[4])
                    except IndexError:
                        chain = self.get_chain(split_line[3])
                    chain.add_atom(split_line)
    def parse_all(self):
        for chain in self.chains:
            chain.parse()
def ramachandran_plot(protein):
    plt.xlabel('$\phi$')
    plt.ylabel('$\psi$')
    plt.axis([-180,180,-180,180])
    major_ticks = np.arange(-180,181,180)
    ax = plt.gca()
    ax.set_xticks(major_ticks)
    ax.set_yticks(major_ticks)
    ax.set_xticks([0],minor = True)
    ax.set_yticks([0],minor = True)
    ax.grid(color = 'grey',linestyle = '-', linewidth = 1, alpha = 0.5)
    plt.title('Ramachandran Plot')
    lines = list()
    labels = list()
    for idx,chain in enumerate(protein.chains):
        phi = [d[0] for d in chain.dihedrals[1:-1]]
        psi = [d[1] for d in chain.dihedrals[1:-1]]
        line = ax.plot(phi,
                       psi,
                       COLOR_PLATTER[idx]+'o',
                       markersize = 4)
        lines.append(line)
        labels.append('Chain %c'%(chr(idx+65)))
    ax.legend(labels)
#    plt.show()

if __name__ == "__main__":
    test_path = "/home/heavens/CMU/Structure_Biology/HW1/test_data/3i40.pdb"
    output_file = "/home/heavens/CMU/Structure_Biology/HW1/HW1_HaotianTeng_HAT65.md"
    output_image = "/home/heavens/CMU/Structure_Biology/HW1/ramachandran.png"
    p = Protein()
    p.read_pdb(test_path)         
    p.parse_all()
    content = mark_down_print_protein(p)
    with open(output_file,'w') as f:
        f.write(content)
        f.write("## Ramachandran Plot  \n")
        f.write("![Ramachandran Plot](%s)"%(os.path.basename(output_image)))
    ramachandran_plot(p)
    plt.savefig(output_image)