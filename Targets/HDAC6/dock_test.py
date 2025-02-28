'''
Author: Rui Qin
Date: 2025-02-28 15:03:19
LastEditTime: 2025-02-28 16:32:24
Description: 

'''
from vina import Vina
import subprocess
import rdkit.Chem as Chem
from rdkit.Chem import AllChem

lig_pdbqt = './5wgl_S_AH4.pdbqt'
rectz_pdbqt = 'HDAC6_8bjk_rec_A_tz.pdbqt'
out_pdbqt = '5wgl_S_AH4_redock.pdbqt'

v = Vina(sf_name='ad4')
#v.set_receptor(rectz_pdbqt)
v.set_ligand_from_file(lig_pdbqt)
v.load_maps('maps/HDAC6_8bjk_rec_A_tz')
v.dock(n_poses=9, exhaustiveness=32)
v.write_pose(out_pdbqt, overwrite=True)