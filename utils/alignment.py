'''
Author: haoqiang haoqiang@mindrank.ai
Date: 2022-07-13 09:08:39
LastEditors: haoqiang haoqiang@mindrank.ai
LastEditTime: 2022-07-29 10:06:24
FilePath: /work-home/molecule-3d-similarity/utils/alignment.py
Description: some tools for aligning molecular conformers

Copyright (c) 2022 by haoqiang haoqiang@mindrank.ai, All Rights Reserved. 
'''

from copy import deepcopy
from joblib import Parallel, delayed

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
def rdkit_o3a(prb_mols: list, ref_mol: Mol, prb_mmff=None, ref_mmff=None, n_jobs=1, verbose=0) -> list:
    
    def o3a_align(prb_mol, ref_mol, prb_mmff, ref_mmff):
        prb_mol = deepcopy(prb_mol)
        o3a = rdMolAlign.GetO3A(prb_mol, ref_mol, prbPyMMFFMolProperties=prb_mmff, refPyMMFFMolProperties=ref_mmff)
        o3a.Align()
        return prb_mol, o3a.Score()

    a_list = Parallel(n_jobs=n_jobs, verbose=verbose)\
                            (delayed(o3a_align)\
                            (prb_mol, ref_mol, prb_mmff[i] if prb_mmff else prb_mmff, ref_mmff,) \
                                for i, prb_mol in enumerate(prb_mols)) 
    
    return [ele[0] for ele in a_list], [ele[1] for ele in a_list]

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
# def rdkit_crippeno3a(query_mol_list: list, ref_mol: Mol, prb_crippen=[], ref_crippen=[]):
#     score_list = []
#     for i, query_mol in enumerate(query_mol_list):
#         o3a = rdMolAlign.GetCrippenO3A(
#             query_mol, 
#             ref_mol,
#             prbCrippenContribs=prb_crippen[i] if prb_crippen else prb_crippen,
#             refCrippenContribs=ref_crippen,
#             )
#         o3a.Align()
#         score_list.append(o3a.Score())
#     return score_list

def rdkit_crippeno3a(prb_mols: list, ref_mol: Mol, prb_crippen=[], ref_crippen=[], n_jobs=1, verbose=0) -> list:
    
    def crippeno3a_align(prb_mol, ref_mol, prb_crippen, ref_crippen):
        prb_mol = deepcopy(prb_mol)
        o3a = rdMolAlign.GetCrippenO3A(prb_mol, ref_mol, prbCrippenContribs=prb_crippen, refCrippenContribs=ref_crippen)
        o3a.Align()
        return prb_mol, o3a.Score()

    a_list = Parallel(n_jobs=n_jobs, verbose=verbose)\
                            (delayed(crippeno3a_align)\
                            (prb_mol, ref_mol, prb_crippen[i] if prb_crippen else prb_crippen, ref_crippen,) \
                                for i, prb_mol in enumerate(prb_mols)) 
    
    return [ele[0] for ele in a_list], [ele[1] for ele in a_list]
    
