'''
Author: haoqiang haoqiang@mindrank.ai
Date: 2022-07-13 14:22:18
LastEditors: haoqiang haoqiang@mindrank.ai
LastEditTime: 2022-08-01 03:57:14
FilePath: /work-home/molecule-3d-similarity/utils/similarity.py
Description: Investigated methods for calculating the similarity of molecular 3D structures

Copyright (c) 2022 by haoqiang haoqiang@mindrank.ai, All Rights Reserved. 
'''
import os
import sys

base_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(base_path)

from joblib import Parallel, delayed
# import dill as pickle

import pandas as pd

from ast import Global
from rdkit import RDConfig
from rdkit.Chem import AllChem, rdShapeHelpers, Get3DDistanceMatrix, TorsionFingerprints, rdMolAlign, MACCSkeys
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate
from rdkit import DataStructs, Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import GetMorganFingerprint
from rdkit.DataStructs import TanimotoSimilarity, TverskySimilarity, FingerprintSimilarity, DiceSimilarity
from rdkit.Chem.FeatMaps import FeatMaps
from oddt.shape import usr, usr_cat, electroshape, usr_similarity
from espsim import GetShapeSim, GetEspSim
from espsim.helpers import mlCharges

from utils.alignment import rdkit_o3a, rdkit_crippeno3a
from utils.functions import calc_SC_score, calc_SC_score_tanimoto, calc_SC_score_tversky, get_esp


'''
description: 
    Compute the shape protrude distance between two molecule based on a predefined alignment
    https://www.rdkit.org/docs/source/rdkit.Chem.rdShapeHelpers.html
param {*} query_mol_list
param {*} ref_mol
param {*} allowReordering
return {*}
'''
def rdkit_shape_protrude_dist(prb_mols, ref_mol, allowReordering=False, n_jobs=1, verbose=0):
    
    def p_func(prb_mol):
        p_dist = rdShapeHelpers.ShapeProtrudeDist(prb_mol, ref_mol, allowReordering=allowReordering)
        
        return 1- p_dist
        
    dist_list = Parallel(n_jobs=n_jobs, verbose=verbose)(delayed(p_func)(prb_mol) for prb_mol in prb_mols)  

    return dist_list


'''
description: 
    Compute the shape tanimoto distance between two molecule based on a predefined alignment
    https://www.rdkit.org/docs/source/rdkit.Chem.rdShapeHelpers.html
param {*} query_mol_list
param {*} ref_mol
return {*}
'''
def rdkit_shape_tanimoto_dist(prb_mols, ref_mol, n_jobs=1, verbose=0):

    def p_func(prb_mol):
        t_dist = rdShapeHelpers.ShapeTanimotoDist(prb_mol, ref_mol)
        return 1- t_dist
    
    dist_list = Parallel(n_jobs=n_jobs, verbose=verbose)(delayed(p_func)(prb_mol) for prb_mol in prb_mols) 

    return dist_list


'''
description: 
    Compute the shape tversky index between two molecule based on a predefined alignment
    https://www.rdkit.org/docs/source/rdkit.Chem.rdShapeHelpers.html
param {*} query_mol_list
param {*} ref_mol
param {*} alpha
param {*} beta
return {*}
'''
def rdkit_shape_tversky_index(prb_mols, ref_mol, alpha=0.7, beta=0.3, n_jobs=1, verbose=0):
    
    def p_func(prb_mol):
        t_index = rdShapeHelpers.ShapeTverskyIndex(prb_mol, ref_mol, alpha=alpha, beta=beta)
        return t_index
    
    index_list = Parallel(n_jobs=n_jobs, verbose=verbose)(delayed(p_func)(prb_mol) for prb_mol in prb_mols) 
    
    return index_list


'''
description: 
    Same Molecules.
    The RDKit provides an implementation of the torsion fingerprint deviation \
    (TFD) approach developed by Schulz-Gasch et al. (J. Chem. Inf. Model, 52, 1499, 2012).
    https://rdkit.readthedocs.io/en/latest/Cookbook.html
param {*} query_mol_list
param {*} ref_mol
return {*}
'''
def rdkit_tfd(query_mol_list, ref_mol):
    tfd_list = []
    for query_mol in query_mol_list:
        tfd = TorsionFingerprints.GetTFDBetweenMolecules(query_mol, ref_mol)
        tfd_list.append(tfd)

    return tfd_list



fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')    
fdef = AllChem.BuildFeatureFactory(fdefName)
fmParams = {}
keep = ('Donor', 'Acceptor', 'NegIonizable', 'PosIonizable',
        'ZnBinder', 'Aromatic', 'Hydrophobe', 'LumpedHydrophobe')
for k in fdef.GetFeatureFamilies():
    fparams = FeatMaps.FeatMapParams()
    fmParams[k] = fparams
def rdkit_sc_score(prb_mols, ref_mol, n_jobs=1, verbose=0):

    def p_func(prb_mol):
        global fdefName
        fdef = AllChem.BuildFeatureFactory(fdefName)
        global keep
        global fmParams

        score = calc_SC_score(
            query_mol=prb_mol,
            ref_mol=ref_mol,
            fdef=fdef,
            keep=keep,
            fmParams=fmParams
        )
        return score

    score_list = Parallel(n_jobs=n_jobs, verbose=verbose)(delayed(p_func)(prb_mol) for prb_mol in prb_mols) 

    return score_list



def rdkit_sc_score_tanimoto(prb_mols, ref_mol, n_jobs=1, verbose=0):
    
    def p_func(prb_mol):
        global fdefName
        fdef = AllChem.BuildFeatureFactory(fdefName)
        global keep
        global fmParams
        
        score = calc_SC_score_tanimoto(
            query_mol=prb_mol,
            ref_mol=ref_mol,
            fdef=fdef,
            keep=keep,
            fmParams=fmParams
        )
        return score

    score_list = Parallel(n_jobs=n_jobs, verbose=verbose)(delayed(p_func)(prb_mol) for prb_mol in prb_mols) 
    
    return score_list


def rdkit_sc_score_tversky(prb_mols, ref_mol, n_jobs=1, verbose=0):

    
    def p_func(prb_mol):
        global fdefName
        fdef = AllChem.BuildFeatureFactory(fdefName)
        global keep
        global fmParams
        
        score = calc_SC_score_tversky(
            query_mol=prb_mol,
            ref_mol=ref_mol,
            fdef=fdef,
            keep=keep,
            fmParams=fmParams
        )
        return score

    score_list = Parallel(n_jobs=n_jobs, verbose=verbose)(delayed(p_func)(prb_mol) for prb_mol in prb_mols) 
    
    return score_list


'''
description: 
    A 3D pharmacophore fingerprint can be calculated using the RDKit by feeding a \
    3D distance matrix to the 2D-pharmacophore machinery.
    https://rdkit.readthedocs.io/en/latest/Cookbook.html
    https://www.rdkit.org/docs/source/rdkit.DataStructs.cDataStructs.html
param {*} query_mol_list
param {*} ref_mol
return {*}
'''
def rdkit_pharm_tanimoto(prb_mols, ref_mol, n_jobs=1, verbose=0):
    
    similarity_list = []
    factory = Gobbi_Pharm2D.factory
    ref_fp = Generate.Gen2DFingerprint(ref_mol, factory, dMat=Chem.Get3DDistanceMatrix(ref_mol))
    
    def p_func(prb_mol):
        factory = Gobbi_Pharm2D.factory
        prb_fp = Generate.Gen2DFingerprint(prb_mol, factory, dMat=Chem.Get3DDistanceMatrix(prb_mol))
        similarity = DataStructs.TanimotoSimilarity(ref_fp, prb_fp)
        return similarity
    
    similarity_list = Parallel(n_jobs=n_jobs, verbose=verbose)(delayed(p_func)(prb_mol) for prb_mol in prb_mols) 

    return similarity_list


'''
description:
    A 3D pharmacophore fingerprint can be calculated using the RDKit by feeding a \
    3D distance matrix to the 2D-pharmacophore machinery.
    https://rdkit.readthedocs.io/en/latest/Cookbook.html
    The Tversky similarity measure is asymmetric..
    https://www.rdkit.org/docs/source/rdkit.DataStructs.cDataStructs.html
param {*} query_mol_list
param {*} ref_mol
param {*} alpha
param {*} beta
return {*}
'''
def rdkit_pharm_tversky(prb_mols, ref_mol, alpha=0.7, beta=0.3, n_jobs=1, verbose=0):
    
    similarity_list = []
    factory = Gobbi_Pharm2D.factory
    ref_fp = Generate.Gen2DFingerprint(ref_mol, factory, dMat=Chem.Get3DDistanceMatrix(ref_mol))
    
    def p_func(prb_mol):
        factory = Gobbi_Pharm2D.factory
        prb_fp = Generate.Gen2DFingerprint(prb_mol, factory, dMat=Chem.Get3DDistanceMatrix(prb_mol))
        similarity = DataStructs.TverskySimilarity(ref_fp, prb_fp, a=alpha, b=beta)
        return similarity
    
    similarity_list = Parallel(n_jobs=n_jobs, verbose=verbose)(delayed(p_func)(prb_mol) for prb_mol in prb_mols) 

    return similarity_list


def rdkit_fp_maccs(prb_mols, ref_mol, n_jobs=1, verbose=0):
    ref_fp = MACCSkeys.GenMACCSKeys(ref_mol)

    def p_func(prb_mol):
        prb_fp = MACCSkeys.GenMACCSKeys(prb_mol)
        similarity = FingerprintSimilarity(ref_fp, prb_fp)
        return similarity

    similarity_list = Parallel(n_jobs=n_jobs, verbose=verbose)(delayed(p_func)(prb_mol) for prb_mol in prb_mols) 

    return similarity_list


def rdkit_fp_maccs_tainimoto(prb_mols, ref_mol, n_jobs=1, verbose=0):
    ref_fp = MACCSkeys.GenMACCSKeys(ref_mol)

    def p_func(prb_mol):
        prb_fp = MACCSkeys.GenMACCSKeys(prb_mol)
        similarity = DataStructs.TanimotoSimilarity(ref_fp, prb_fp)
        return similarity

    similarity_list = Parallel(n_jobs=n_jobs, verbose=verbose)(delayed(p_func)(prb_mol) for prb_mol in prb_mols) 

    return similarity_list


def rdkit_fp_maccs_tversky(prb_mols, ref_mol, alpha=0.7, beta=0.3, n_jobs=1, verbose=0):
    ref_fp = MACCSkeys.GenMACCSKeys(ref_mol)

    def p_func(prb_mol):
        prb_fp = MACCSkeys.GenMACCSKeys(prb_mol)
        similarity = DataStructs.TverskySimilarity(ref_fp, prb_fp, a=alpha, b=beta)
        return similarity

    similarity_list = Parallel(n_jobs=n_jobs, verbose=verbose)(delayed(p_func)(prb_mol) for prb_mol in prb_mols) 

    return similarity_list


def rdkit_fp_morgan(prb_mols, ref_mol, radius=2, n_jobs=1, verbose=0):
    ref_fp = AllChem.GetMorganFingerprint(ref_mol, radius)

    def p_func(prb_mol):
        prb_fp = AllChem.GetMorganFingerprint(prb_mol, radius)
        similarity = DataStructs.DiceSimilarity(ref_fp, prb_fp)
        return similarity

    similarity_list = Parallel(n_jobs=n_jobs, verbose=verbose)(delayed(p_func)(prb_mol) for prb_mol in prb_mols) 

    return similarity_list


def rdkit_fp_morgan_tanimoto(prb_mols, ref_mol, radius=2, n_jobs=1, verbose=0):
    ref_fp = AllChem.GetMorganFingerprint(ref_mol, radius)

    def p_func(prb_mol):
        prb_fp = AllChem.GetMorganFingerprint(prb_mol, radius)
        similarity = DataStructs.TanimotoSimilarity(ref_fp, prb_fp)
        return similarity

    similarity_list = Parallel(n_jobs=n_jobs, verbose=verbose)(delayed(p_func)(prb_mol) for prb_mol in prb_mols) 

    return similarity_list

def rdkit_fp_morgan_tversky(prb_mols, ref_mol, radius=2, alpha=0.7, beta=0.3, n_jobs=1, verbose=0):
    ref_fp = AllChem.GetMorganFingerprint(ref_mol, radius)

    def p_func(prb_mol):
        prb_fp = AllChem.GetMorganFingerprint(prb_mol, radius)
        similarity = DataStructs.TverskySimilarity(ref_fp, prb_fp, a=alpha, b=beta)
        return similarity

    similarity_list = Parallel(n_jobs=n_jobs, verbose=verbose)(delayed(p_func)(prb_mol) for prb_mol in prb_mols) 

    return similarity_list


'''
description: 
    Ballester PJ, Richards WG (2007). Ultrafast shape recognition to search \
    compound databases for similar molecular shapes. Journal of computational \
    chemistry, 28(10):1711-23. 
    http://dx.doi.org/10.1002/jcc.20681
    https://oddt.readthedocs.io/en/latest/index.html#usage-instructions
param {*} query_mol_list
param {*} ref_mol
return {*}
'''
def oddt_usr(prb_mols, ref_mol, n_jobs=1, verbose=0):
    ref_shape = usr(ref_mol)

    def p_func(prb_mol):
        prb_shape = usr(prb_mol)
        similarity = usr_similarity(ref_shape, prb_shape)
        return similarity
    
    similarity_list = Parallel(n_jobs=n_jobs, verbose=verbose)(delayed(p_func)(prb_mol) for prb_mol in prb_mols) 
    
    return similarity_list


'''
description: 
    USRCAT: real-time ultrafast shape recognition with pharmacophoric constraints
    http://dx.doi.org/10.1186/1758-2946-4-27
    https://oddt.readthedocs.io/en/latest/index.html#usage-instructions
param {*} query_mol_list
param {*} ref_mol
return {*}
'''
def oddt_usr_cat(prb_mols, ref_mol, n_jobs=1, verbose=0):
    ref_shape = usr_cat(ref_mol)

    def p_func(prb_mol):
        prb_shape = usr_cat(prb_mol)
        similarity = usr_similarity(ref_shape, prb_shape)
        return similarity
    
    similarity_list = Parallel(n_jobs=n_jobs, verbose=verbose)(delayed(p_func)(prb_mol) for prb_mol in prb_mols) 
    
    return similarity_list

'''
description: 
    ElectroShape: fast molecular similarity calculations incorporating \
    shape, chirality and electrostatics
    http://dx.doi.org/doi:10.1007/s10822-010-9374-0
    https://oddt.readthedocs.io/en/latest/index.html#usage-instructions
param {*} query_mol_list
param {*} ref_mol
return {*}
'''
def oddt_electroshape(prb_mols, ref_mol, n_jobs=1, verbose=0):
    ref_shape = electroshape(ref_mol)

    def p_func(prb_mol):
        prb_shape = electroshape(prb_mol)
        similarity = usr_similarity(ref_shape, prb_shape)
        return similarity
    
    similarity_list = Parallel(n_jobs=n_jobs, verbose=verbose)(delayed(p_func)(prb_mol) for prb_mol in prb_mols) 
    
    return similarity_list


def acpc(ref_path, prb_path, out_path, num_core):
    assert (ref_path.split('.')[-1] == prb_path.split('.')[-1]) or (prb_path.split('.')[-1]=='bin'), 'make sure ref and prb file have same format'
    if not os.path.exists(os.path.split(out_path)[0]):
        os.system(f'mkdir -p {os.path.split(out_path)[0]}')
    os.system(f"acpc_par -np {num_core} -q {ref_path} -db {prb_path} -o {out_path}")

    out_df = pd.read_csv(out_path, sep=' ', header=None, names=['name', 'id', 'score'])
    out_df['convert'] = out_df.apply(lambda x: (x['score'], 1, x['id']) if 'actives' in x['name'] else (x['score'], 0, x['id']), axis=1)

    return out_df['convert'].to_list()


def esp_sim(prb_mols, ref_mol, align='done', shape_sim='tanimoto', systems='mmff', metric='tanimoto', prb_charge=None, ref_charge=None, n_jobs=1, verbose=0):
    # align molecules 
    if   align == 'o3a':
        prb_mols, _ = rdkit_o3a(prb_mols=prb_mols, ref_mol=ref_mol, n_jobs=n_jobs, verbose=verbose)
    elif align == 'crippeno3a':
        prb_mols, _ = rdkit_crippeno3a(prb_mols=prb_mols, ref_mol=ref_mol, n_jobs=n_jobs, verbose=verbose)

    def p_func(prb_mol, prb_charge=None):

        # calculate shape and esp similarities
        if   shape_sim == 'protrude':
            shape = 1 - AllChem.ShapeProtrudeDist(prb_mol, ref_mol)
        elif shape_sim == 'tanimoto':
            shape = 1 - AllChem.ShapeTanimotoDist(prb_mol, ref_mol)
        elif shape_sim == 'tversky':
            shape = AllChem.ShapeTverskyIndex(prb_mol, ref_mol, alpha=0.7, beta=0.3)
        else:
            shape = 0

        if   systems == 'mmff':
            esp = get_esp(prb_mol, ref_mol, 0, 0, 'mmff',      metric=metric)
        elif systems == 'gasteiger':
            esp = get_esp(prb_mol, ref_mol, 0, 0, 'gasteiger', metric=metric)
        elif systems == 'ml':
            esp = get_esp(prb_mol, ref_mol, 0, 0, 'ml',        metric=metric, prbCharge=prb_charge, refCharge=ref_charge)
            
        return shape+esp

    similarity_list = Parallel(n_jobs=n_jobs, verbose=verbose)(delayed(p_func)(prb_mol, prb_charge[i] if prb_charge else None) for i, prb_mol in enumerate(prb_mols)) 
    
    return similarity_list