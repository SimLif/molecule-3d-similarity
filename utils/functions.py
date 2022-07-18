'''
Author: haoqiang haoqiang@mindrank.ai
Date: 2022-07-13 09:53:29
LastEditors: haoqiang haoqiang@mindrank.ai
LastEditTime: 2022-07-15 10:35:46
FilePath: /work-home/molecule-3d-similarity/utils/functions.py
Description: some useful toolkits for molecule

Copyright (c) 2022 by haoqiang haoqiang@mindrank.ai, All Rights Reserved. 
'''
import os
import time

from rdkit import Chem
import nglview as nv
from openbabel import openbabel
from rdkit import RDConfig
from rdkit.Chem import AllChem, rdShapeHelpers, Get3DDistanceMatrix, TorsionFingerprints, rdMolAlign
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate
from rdkit.DataStructs import TanimotoSimilarity, TverskySimilarity
from rdkit.Chem.FeatMaps import FeatMaps


'''
description: covert molecular formation
param {*} from_file_path
param {*} from_format
param {*} to_format
return {*}
'''
def convert_format(from_file_path, from_format, to_format):
    # define
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats(from_format, to_format)

    # read
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, from_file_path)

    # correct
    mol.AddHydrogens()

    # save
    to_file_path = from_file_path.split('.')[0] + '.' + to_format
    obConversion.WriteFile(mol, to_file_path)

    return to_file_path


'''
description: read mol2 file which contains multiple molecules
param {*} file
param {*} sanitize
return {*}
'''
def Mol2MolSupplier(file=None ,sanitize=True, removeHs=True):
    mols=[]
    with open(file, 'r') as f:
        doc=[line for line in f.readlines()]

    mark=[index for (index,p) in enumerate(doc) if '@<TRIPOS>MOLECULE' in p]
    mark.append(len(doc)+1)
    
    interval=list(zip(mark[:-1],mark[1:]))
    for i in interval:
        block = "".join(doc[i[0]:i[1]])
        m=Chem.MolFromMol2Block(block, sanitize=sanitize, removeHs=removeHs)
        if m:
            mols.append(m)
    return mols


"""Generate a view of the ligand molecules.

Parameters
-----------
molecules: list of rdkit.Chem.rdchem.Mol

Returns
----------
nglview.widget.NGLWidget
"""
def show_ligands(molecules):
    view = nv.NGLWidget()
    component = view.add_component(molecules[0])
    time.sleep(0.1)
    component.clear()
    # component.add_ball_and_stick(multipleBond=True, color='green')
    component.add_surface(color='green', opacity=0.3)
    for molecule in molecules[1:]:
        component = view.add_component(molecule)
        time.sleep(0.1)
        component.clear()
        component.add_surface(opacity=0.3)
    return view


'''
description: 
param {*} query_mol
param {*} ref_mol
param {*} fdef
param {*} keep
param {*} fmParams
return {*}
'''
def get_FeatureMapScore(query_mol, ref_mol, fdef, keep, fmParams):
    featLists = []
    for m in [query_mol, ref_mol]:
        rawFeats = fdef.GetFeaturesForMol(m)
        # filter that list down to only include the ones we're intereted in
        featLists.append([f for f in rawFeats if f.GetFamily() in keep])
    fms = [FeatMaps.FeatMap(feats=x, weights=[1] * len(x), params=fmParams) for x in featLists]
    fms[0].scoreMode = FeatMaps.FeatMapScoreMode.Best
    fm_score = fms[0].ScoreFeats(featLists[1]) / min(fms[0].GetNumFeatures(), len(featLists[1]))

    return fm_score


'''
description: 
param {*} query_mol
param {*} ref_mol
param {*} fdef
param {*} keep
param {*} fmParams
return {*}
'''
def calc_SC_score(query_mol, ref_mol, fdef, keep, fmParams):
    # _ = rdMolAlign.GetO3A(query_mol, ref_mol).Align()
    fm_score = get_FeatureMapScore(query_mol, ref_mol, fdef, keep, fmParams)

    protrude_dist = rdShapeHelpers.ShapeProtrudeDist(query_mol, ref_mol,
                                                     allowReordering=False)
    SC_RDKit_score = 0.5 * fm_score + 0.5 * (1 - protrude_dist)

    return SC_RDKit_score


'''
description: 
    use ShapeTanimotoDist
param {*} query_mol
param {*} ref_mol
param {*} fdef
param {*} keep
param {*} fmParams
return {*}
'''
def calc_SC_score_plus(query_mol, ref_mol, fdef, keep, fmParams):
    # _ = rdMolAlign.GetCrippenO3A(query_mol, ref_mol).Align()
    fm_score = get_FeatureMapScore(query_mol, ref_mol, fdef, keep, fmParams)

    protrude_dist = rdShapeHelpers.ShapeTanimotoDist(query_mol, ref_mol)
    SC_RDKit_score = 0.5 * fm_score + 0.5 * (1 - protrude_dist)

    return SC_RDKit_score


'''
description: 
    use TverskyIndex
param {*} query_mol
param {*} ref_mol
param {*} fdef
param {*} keep
param {*} fmParams
return {*}
'''
def calc_SC_score_tversky(query_mol, ref_mol, fdef, keep, fmParams, alpha=0.7, beta=0.3):
    # _ = rdMolAlign.GetCrippenO3A(query_mol, ref_mol).Align()
    fm_score = get_FeatureMapScore(query_mol, ref_mol, fdef, keep, fmParams)

    tversky_index = rdShapeHelpers.ShapeTverskyIndex(
            query_mol, 
            ref_mol,
            alpha=alpha,
            beta=beta,
            )
    SC_RDKit_score = 0.5 * fm_score + 0.5 * tversky_index

    return SC_RDKit_score
