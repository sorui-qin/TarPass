'''
Author: Rui Qin
Date: 2025-03-01 14:58:38
LastEditTime: 2025-03-05 21:01:02
Description: 
'''
from vina import Vina
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy
from meeko import PDBQTMolecule
from meeko import RDKitMolCreate
from rdkit import Chem
import glob
import tqdm
import json
from pathlib import Path
import numpy as np


def sdf2centroid(sdf_file):
    supp = Chem.SDMolSupplier(sdf_file, sanitize=False)
    lig_xyz = supp[0].GetConformer().GetPositions()
    centroid_x = lig_xyz[:,0].mean()
    centroid_y = lig_xyz[:,1].mean()
    centroid_z = lig_xyz[:,2].mean()
    return [centroid_x, centroid_y, centroid_z]


def lig_prep(lig_sdf):
    mol = Chem.SDMolSupplier(lig_sdf, removeHs=False)[0]
    mk_prep = MoleculePreparation()
    molsetup_list = mk_prep(mol)
    molsetup = molsetup_list[0]
    pdbqt_string = PDBQTWriterLegacy.write_string(molsetup)[0]
    return pdbqt_string


def maps_prep(pocket, box_size=[20, 20, 20]):
    assert pocket != 'HDAC6', 'Preparation for HDAC6 should follow `Target/HDAC6/preprocess.ipynb`'

    rec_pdbqt = glob.glob(f'dock/pdbqtfiles/{pocket}_*.pdbqt')[0]
    try: # for apo structure
        ref_sdf = glob.glob(f'Targets/{pocket}/*ligand*.sdf')[0]
        center = sdf2centroid(ref_sdf)
    except: # for holo and AlphaFold structure
        center_json = glob.glob(f'Targets/{pocket}/center.json')[0]
        with open(center_json, "r") as f:
            center = json.load(f)["center"]

    Path(f'dock/maps/{pocket}').mkdir(exist_ok=True)

    v = Vina(cpu=0, seed=42)
    v.set_receptor(rec_pdbqt)
    v.compute_vina_maps(center, box_size)
    v.write_maps(map_prefix_filename=f'dock/maps/{pocket}/{pocket}', overwrite=True)


def dock(pdbqt_string, maps, mode='dock', sf_name='vina', seed=0, exhaust=32, n_poses=1, verbose=0) -> tuple[float, list]:
    """_summary_

    Args:
        pdbqt_string (_type_): _description_
        maps (_type_): _description_
        mode (str, optional): _description_. Defaults to 'dock'.
        sf_name (str, optional): _description_. Defaults to 'vina'.
        seed (int, optional): _description_. Defaults to 0.
        exhaust (int, optional): _description_. Defaults to 32.
        n_poses (int, optional): _description_. Defaults to 1.

    Returns:
        tuple[float, list]: _description_
    """
    v = Vina(sf_name=sf_name, cpu=0, seed=seed, verbosity=verbose)
    v.set_ligand_from_string(pdbqt_string)
    v.load_maps(maps)
    if mode == 'score_only':
        score = v.score()[0]
        return score, None
    if mode == 'dock':
        v.dock(exhaustiveness=exhaust, n_poses=n_poses)
        score = v.energies(n_poses=n_poses)[0][0]
        vina_output_string = v.poses(n_poses=n_poses)
        pdbqt_mol = PDBQTMolecule(vina_output_string, skip_typing=True)
        poses = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
        return score, poses
    

def write_poses(poses, filename):
    with Chem.SDWriter(filename) as w:
        for mol in poses:
            w.write(mol)


def cal_RMSD(ref_mol, check_mol):
    rmsd = Chem.rdMolAlign.CalcRMS(
        Chem.RemoveHs(check_mol), Chem.RemoveHs(ref_mol))
    return rmsd


if __name__ == '__main__':
    target_dirs = sorted(glob.glob('./Targets/*/'))
    target_dirs.remove('./Targets/5HT2A-AF/')
    target_dirs.remove('./Targets/BRD4-holo/')

    for target in tqdm.tqdm(target_dirs):
        pocket = Path(target).name

        #ligprep
        ref_sdf = glob.glob(f'{target}/*ligand*.sdf')[0]
        pdbqt_string = lig_prep(ref_sdf)

        #recprep
        maps_dir = f'dock/maps/{pocket}/'
        if not Path(maps_dir).exists():
            maps_prep(pocket)
        
        #dock
        assert Path(maps_dir).exists(), 'Affinity map files do not exist!'
        outfile = f'{target}/redock-32.sdf'
        if not Path(outfile).exists():
            sf_name = 'ad4' if pocket == 'HDAC6' else 'vina'
            maps = f'dock/maps/{pocket}/{pocket}'
            _, poses = dock(pdbqt_string=pdbqt_string,
                            maps=maps,
                            sf_name=sf_name,
                            seed=42,
                            verbose=1)
            with Chem.SDWriter(outfile) as w:
                for mol in poses:
                    w.write(mol)