'''
Author: Rui Qin
Date: 2025-03-01 14:58:38
LastEditTime: 2025-03-05 09:20:47
Description: 
'''
from vina import Vina
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy
from meeko import PDBQTMolecule
from meeko import RDKitMolCreate
from rdkit import Chem
import subprocess
import glob
import tqdm
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


def rec_prep(rec_prefix):
    command_prep = f'''mk_prepare_receptor.py \
                    -i {rec_prefix}.pdb \
                    -o {rec_prefix} -p'''
    proc = subprocess.Popen(command_prep, shell=True, stdin=subprocess.PIPE, 
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout_data, stderr_data = proc.communicate()
    print(stdout_data.decode(), stderr_data.decode())


def dock(pdbqt_string, rec_prefix, center):
    box_size = [20, 20, 20]
    v = Vina(cpu=0, seed=42)
    v.set_receptor(f'{rec_prefix}.pdbqt')
    v.set_ligand_from_string(pdbqt_string)
    v.compute_vina_maps(center, box_size)
    v.dock(exhaustiveness=32, n_poses=10)
    vina_output_string = v.poses(n_poses=10)
    pdbqt_mol = PDBQTMolecule(vina_output_string, skip_typing=True)
    rdkitmol_list = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
    return rdkitmol_list


def dock_zinc(pdbqt_string, maps):
    v = Vina(sf_name='ad4',seed=42)
    v.set_ligand_from_string(pdbqt_string)
    v.load_maps(maps)
    v.dock(exhaustiveness=32, n_poses=10)
    vina_output_string = v.poses(n_poses=10)
    pdbqt_mol = PDBQTMolecule(vina_output_string, skip_typing=True)
    rdkitmol_list = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
    return rdkitmol_list


def cal_RMSD(ref_mol, check_mol):
    rmsd = Chem.rdMolAlign.CalcRMS(
        Chem.RemoveHs(check_mol), Chem.RemoveHs(ref_mol))
    return rmsd


if __name__ == '__main__':
    target_dirs = sorted(glob.glob('./Targets/*/'))
    target_dirs.remove('./Targets/5HT2A_AF/')
    target_dirs.remove('./Targets/BRD4_holo/')
    target_dirs.remove('./Targets/HDAC6/')

    rmsd_list = []
    for target in tqdm.tqdm(target_dirs):
        ref_sdf = glob.glob(f'{target}/*ligand*.sdf')[0]
        rec_pdbs = glob.glob(f'{target}/*rec*.pdb')
        if len(rec_pdbs) == 1:
            rec_pdb = rec_pdbs[0]
        else:
            rec_pdb = [rec for rec in rec_pdbs 
                       if not rec.endswith('empty.pdb')][0]
        rec_prefix = rec_pdb.replace('.pdb', '')
        redock_sdf = f'{target}/redock_32-multi.sdf'

        # Ligprep
        pdbqt_string = lig_prep(ref_sdf)

        # Recprep
        if not Path(f'{rec_prefix}.pdbqt').exists():
            rec_prep(rec_prefix)
        
        # Docking
        assert Path(f'{rec_prefix}.pdbqt').exists()
        if not Path(redock_sdf).exists():
            center = sdf2centroid(ref_sdf)
            docked_mols  = dock(pdbqt_string, rec_prefix, center)
        
            # Calculate RMSD
            #ref_mol = Chem.SDMolSupplier(ref_sdf)[0]
            #rmsd_list.append(cal_RMSD(ref_mol, docked_mol))
            with Chem.SDWriter(redock_sdf) as w:
                for mol in docked_mols:
                    w.write(mol)
    
    # Zinc dock
    target = './Targets/HDAC6/'
    ref_sdf = glob.glob(f'{target}/*ligand*.sdf')[0]
    
    pdbqt_string = lig_prep(ref_sdf)
    rec_pdb = glob.glob(f'{target}/*rec*.pdb')[0]

    rec_prefix = rec_pdb.replace('.pdb', '')
    dock_mols = dock_zinc(pdbqt_string, maps=f'{target}/maps/HDAC6_8bjk_rec_A_tz')

    redock_sdf = f'{target}/redock_32-multi.sdf'
    with Chem.SDWriter(redock_sdf) as w:
        for mol in dock_mols:
            w.write(mol)