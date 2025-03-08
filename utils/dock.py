'''
Author: Rui Qin
Date: 2025-03-07 19:49:34
LastEditTime: 2025-03-07 21:05:50
Description: 
'''
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy


def sdf2centroid(sdf_file):
    supp = Chem.SDMolSupplier(sdf_file, sanitize=False)
    lig_xyz = supp[0].GetConformer().GetPositions()
    centroid_x = lig_xyz[:,0].mean()
    centroid_y = lig_xyz[:,1].mean()
    centroid_z = lig_xyz[:,2].mean()
    return [centroid_x, centroid_y, centroid_z]


class LigPrep():
    """Preparation of ligands to PDBQT file
    
    Args:
        input_seq (str): Filename of a sdf file (Path) or a SMILES sequenece.
        optimize (bool): Optimize the conformation (default: 0; Do not optimize)
    """
    def __init__(self, input_seq, optimize=0):
        if input_seq.endswith('.sdf'):
            self.mol = Chem.SDMolSupplier(input_seq, removeHs=False)[0]
            self.optimize = optimize
        else:
            mol = Chem.MolFromSmiles(input_seq)
            if mol:
                self.mol = mol
                self.optimize = True
            else:
                raise ValueError(f'Invalid input, only SDF file or SMILES seq permitted.')

    def ligprep(self):
        self.mol = Chem.AddHs(self.mol, addCoords=True)
        if self.optimize:
            AllChem.EmbedMolecule(self.mol, AllChem.ETKDGv3())
            AllChem.MMFFOptimizeMolecule(self.mol)
        AllChem.ComputeGasteigerCharges(self.mol)

    def get_pdbqt(self):
        self.ligprep()
        mk_prep = MoleculePreparation()
        molsetup_list = mk_prep(self.mol)
        molsetup = molsetup_list[0]
        pdbqt_string = PDBQTWriterLegacy.write_string(molsetup)[0]
        return pdbqt_string