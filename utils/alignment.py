'''
Author: haoqiang haoqiang@mindrank.ai
Date: 2022-07-13 09:08:39
LastEditors: haoqiang haoqiang@mindrank.ai
LastEditTime: 2022-07-13 12:16:24
FilePath: /work-home/molecule-3d-similarity/utils/alignment.py
Description: some tools for aligning molecular conformers

Copyright (c) 2022 by haoqiang haoqiang@mindrank.ai, All Rights Reserved. 
'''

from rdkit.Chem import rdMolAlign, AllChem
from rdkit.Chem.AllChem import Mol


'''
description: 
    Optimally (minimum RMSD) align a molecule to another molecule.
    https://www.rdkit.org/docs/source/rdkit.Chem.rdMolAlign.html
param {list} query_mol_list
param {Mol} ref_mol
return {list} rmsd_list
'''
def rdkit_alignmol(query_mol_list: list, ref_mol: Mol) -> list:
    rmsd_list = []
    for query_mol in query_mol_list:
        rmsd = rdMolAlign.AlignMol(query_mol, ref_mol)
        rmsd_list.append(rmsd)
    return rmsd_list


'''
description: 
    Get an O3A object with atomMap and weights vectors to overlay \
    the probe molecule onto the reference molecule based on MMFF \
    atom types and charges
    https://www.rdkit.org/docs/source/rdkit.Chem.rdMolAlign.html
param {list} query_mol_list
param {Mol} ref_mol
param {*} prb_mmff
param {*} ref_mmff
return {list} score_list
'''
def rdkit_o3a(query_mol_list: list, ref_mol: Mol, prb_mmff=None, ref_mmff=None) -> list:
    score_list = []
    for i, query_mol in enumerate(query_mol_list):
        o3a = rdMolAlign.GetO3A(
            query_mol, 
            ref_mol,
            prbPyMMFFMolProperties=prb_mmff[i] if prb_mmff else prb_mmff,
            refPyMMFFMolProperties=ref_mmff,
            )
        o3a.Align()
        score_list.append(o3a.Score())
    return score_list


'''
description: 
    Get an O3A object with atomMap and weights vectors to overlay \
    the probe molecule onto the reference molecule based on Crippen \
    logP atom contributions
    https://www.rdkit.org/docs/source/rdkit.Chem.rdMolAlign.html
param {list} query_mol_list
param {Mol} ref_mol
param {*} prb_crippen
param {*} ref_crippen
return {list} score list
'''
def rdkit_crippeno3a(query_mol_list: list, ref_mol: Mol, prb_crippen=[], ref_crippen=[]):
    score_list = []
    for i, query_mol in enumerate(query_mol_list):
        o3a = rdMolAlign.GetCrippenO3A(
            query_mol, 
            ref_mol,
            prbCrippenContribs=prb_crippen[i] if prb_crippen else prb_crippen,
            refCrippenContribs=ref_crippen,
            )
        o3a.Align()
        score_list.append(o3a.Score())
    return score_list
    
