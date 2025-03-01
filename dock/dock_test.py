'''
Author: Rui Qin
Date: 2025-03-01 14:58:38
LastEditTime: 2025-03-01 20:13:32
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
    print(stdout_data, stderr_data)


def dock(center, pdbqt_string, rec_prefix):
    box_size = [20, 20, 20]
    v = Vina(cpu=0)
    v.set_receptor(f'{rec_prefix}.pdbqt')
    v.set_ligand_from_string(pdbqt_string)
    v.compute_vina_maps(center, box_size)
    v.dock(exhaustiveness=8, n_poses=1)
    vina_output_string = v.poses(n_poses=1)
    pdbqt_mol = PDBQTMolecule(vina_output_string, skip_typing=True)
    rdkitmol_list = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
    return rdkitmol_list[0]

def cal_RMSD(ref_mol, check_mol):
    rmsd = Chem.rdMolAlign.CalcRMS(
        Chem.RemoveHs(check_mol), Chem.RemoveHs(ref_mol))
    return rmsd

#v = Vina(sf_name='ad4')
#v.set_receptor(prot_pdbqt)
#v.set_ligand_from_file(lig_pdbqt)
#v.load_maps(maps)
#print(v.score()[0])
#v.dock(exhaustiveness=8, n_poses=10)
#v.write_poses('Targets/HDAC6/redock.pdbqt', overwrite=True)

if __name__ == '__main__':
    target_dirs = sorted(glob.glob('./Targets/*/'))
    target_dirs.remove('./Targets/5HT2A_AF/')
    target_dirs.remove('./Targets/BRD4_holo/')
    target_dirs.remove('./Targets/HDAC6/')
    target_dirs.remove('./Targets/HIV-RT/')

    rmsd_list = []
    for target in tqdm.tqdm(target_dirs):
        ref_sdf = glob.glob(f'{target}/*ligand*.sdf')[0]
        rec_pdb = glob.glob(f'{target}/*rec*.pdb')[0]
        rec_prefix = rec_pdb.replace('.pdb', '')

        # Ligprep
        pdbqt_string = lig_prep(ref_sdf)

        # Recprep
        if not Path(f'{rec_prefix}.pdbqt').exists():
            rec_prep(rec_prefix)
        
        # Docking
        assert Path(f'{rec_prefix}.pdbqt').exists()
        center = sdf2centroid(ref_sdf)
        docked_mol = dock(center, pdbqt_string, rec_prefix)
        
        # Calculate RMSD
        ref_mol = Chem.SDMolSupplier(ref_sdf)[0]
        rmsd_list.append(cal_RMSD(ref_mol, docked_mol))
    
    print(f'RMSD MAX: {np.max(rmsd_list)}')
    print(f'RMSD MIN: {np.min(rmsd_list)}')
    print(f'RMSD AVG: {np.mean(rmsd_list)}')
    print(f'RMSD MED: {np.median(rmsd_list)}')