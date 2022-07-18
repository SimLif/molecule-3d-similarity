'''
Author: haoqiang haoqiang@mindrank.ai
Date: 2022-07-13 14:22:18
LastEditors: haoqiang haoqiang@mindrank.ai
LastEditTime: 2022-07-15 10:36:26
FilePath: /work-home/molecule-3d-similarity/utils/similarity.py
Description: Investigated methods for calculating the similarity of molecular 3D structures

Copyright (c) 2022 by haoqiang haoqiang@mindrank.ai, All Rights Reserved. 
'''
import os
import sys

base_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(base_path)

from ast import Global
from rdkit import RDConfig
from rdkit.Chem import AllChem, rdShapeHelpers, Get3DDistanceMatrix, TorsionFingerprints, rdMolAlign
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate
from rdkit.DataStructs import TanimotoSimilarity, TverskySimilarity
from rdkit.Chem.FeatMaps import FeatMaps
from oddt.shape import usr, usr_cat, electroshape, usr_similarity

from utils.functions import calc_SC_score, calc_SC_score_plus, calc_SC_score_tversky


'''
description: 
    Compute the shape protrude distance between two molecule based on a predefined alignment
    https://www.rdkit.org/docs/source/rdkit.Chem.rdShapeHelpers.html
param {*} query_mol_list
param {*} ref_mol
param {*} allowReordering
return {*}
'''
def rdkit_shape_protrude_dist(query_mol_list, ref_mol, allowReordering=False):
    dist_list = []
    for query_mol in query_mol_list:
        protrude_dist = rdShapeHelpers.ShapeProtrudeDist(
            query_mol, 
            ref_mol,
            allowReordering=allowReordering
            )
        dist_list.append(1-protrude_dist)

    return dist_list


'''
description: 
    Compute the shape tanimoto distance between two molecule based on a predefined alignment
    https://www.rdkit.org/docs/source/rdkit.Chem.rdShapeHelpers.html
param {*} query_mol_list
param {*} ref_mol
return {*}
'''
def rdkit_shape_tanimoto_dist(query_mol_list, ref_mol):
    dist_list = []
    for query_mol in query_mol_list:
        tanimoto_dist = rdShapeHelpers.ShapeTanimotoDist(
            query_mol, 
            ref_mol,
            )
        dist_list.append(1-tanimoto_dist)

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
def rdkit_shape_tversky_index(query_mol_list, ref_mol, alpha=0.7, beta=0.3):
    index_list = []
    for query_mol in query_mol_list:
        tversky_index = rdShapeHelpers.ShapeTverskyIndex(
            query_mol, 
            ref_mol,
            alpha=alpha,
            beta=beta,
            )
        index_list.append(tversky_index)

    return index_list


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
def rdkit_pharm_tanimoto(query_mol_list, ref_mol):
    similarity_list = []
    factory = Gobbi_Pharm2D.factory
    ref_fp = Generate.Gen2DFingerprint(ref_mol, factory, dMat=Get3DDistanceMatrix(ref_mol))
    for query_mol in query_mol_list:
        query_fp = Generate.Gen2DFingerprint(query_mol, factory, dMat=Get3DDistanceMatrix(query_mol))
        similarity = TanimotoSimilarity(ref_fp, query_fp)
        similarity_list.append(similarity)

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
def rdkit_pharm_tversky(query_mol_list, ref_mol, alpha=0.7, beta=0.3):
    similarity_list = []
    factory = Gobbi_Pharm2D.factory
    ref_fp = Generate.Gen2DFingerprint(ref_mol, factory, dMat=Get3DDistanceMatrix(ref_mol))
    for query_mol in query_mol_list:
        query_fp = Generate.Gen2DFingerprint(query_mol, factory, dMat=Get3DDistanceMatrix(query_mol))
        similarity = TverskySimilarity(ref_fp, query_fp, a=alpha, b=beta)
        similarity_list.append(similarity)

    return similarity_list


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
def oddt_usr(query_mol_list, ref_mol):
    similarity_list = []
    ref_shape = usr(ref_mol)
    for query_mol in query_mol_list:
        query_shape = usr(query_mol)
        similarity = usr_similarity(ref_shape, query_shape)
        similarity_list.append(similarity)
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
def oddt_usr_cat(query_mol_list, ref_mol):
    similarity_list = []
    ref_shape = usr_cat(ref_mol)
    for query_mol in query_mol_list:
        query_shape = usr_cat(query_mol)
        similarity = usr_similarity(ref_shape, query_shape)
        similarity_list.append(similarity)
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
def oddt_electroshape(query_mol_list, ref_mol):
    similarity_list = []
    ref_shape = electroshape(ref_mol)
    for query_mol in query_mol_list:
        query_shape = electroshape(query_mol)
        similarity = usr_similarity(ref_shape, query_shape)
        similarity_list.append(similarity)
    return similarity_list




def rdkit_sc_score(query_mol_list, ref_mol):
    fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')    
    fdef = AllChem.BuildFeatureFactory(fdefName)
    fmParams = {}
    keep = ('Donor', 'Acceptor', 'NegIonizable', 'PosIonizable',
            'ZnBinder', 'Aromatic', 'Hydrophobe', 'LumpedHydrophobe')
    for k in fdef.GetFeatureFamilies():
        fparams = FeatMaps.FeatMapParams()
        fmParams[k] = fparams
    
    score_list = []
    for query_mol in query_mol_list:
        score = calc_SC_score(
            query_mol=query_mol,
            ref_mol=ref_mol,
            fdef=fdef,
            keep=keep,
            fmParams=fmParams
        )
        score_list.append(score)
    return score_list



def rdkit_sc_score_plus(query_mol_list, ref_mol):
    fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')    
    fdef = AllChem.BuildFeatureFactory(fdefName)
    fmParams = {}
    keep = ('Donor', 'Acceptor', 'NegIonizable', 'PosIonizable',
            'ZnBinder', 'Aromatic', 'Hydrophobe', 'LumpedHydrophobe')
    for k in fdef.GetFeatureFamilies():
        fparams = FeatMaps.FeatMapParams()
        fmParams[k] = fparams
    
    score_list = []
    for query_mol in query_mol_list:
        score = calc_SC_score_plus(
            query_mol=query_mol,
            ref_mol=ref_mol,
            fdef=fdef,
            keep=keep,
            fmParams=fmParams
        )
        score_list.append(score)
    return score_list


def rdkit_sc_score_tversky(query_mol_list, ref_mol):
    fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')    
    fdef = AllChem.BuildFeatureFactory(fdefName)
    fmParams = {}
    keep = ('Donor', 'Acceptor', 'NegIonizable', 'PosIonizable',
            'ZnBinder', 'Aromatic', 'Hydrophobe', 'LumpedHydrophobe')
    for k in fdef.GetFeatureFamilies():
        fparams = FeatMaps.FeatMapParams()
        fmParams[k] = fparams
    
    score_list = []
    for query_mol in query_mol_list:
        score = calc_SC_score_tversky(
            query_mol=query_mol,
            ref_mol=ref_mol,
            fdef=fdef,
            keep=keep,
            fmParams=fmParams
        )
        score_list.append(score)
    return score_list