'''
Author: Rui Qin
Date: 2025-03-07 19:49:34
LastEditTime: 2025-03-13 20:17:58
Description: 
'''
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import pybel
from meeko import MoleculePreparation, PDBQTWriterLegacy
from utils.io import standard_mol, write_sdf


def sdf2centroid(sdf_file):
    supp = Chem.SDMolSupplier(sdf_file, sanitize=False)
    lig_xyz = supp[0].GetConformer().GetPositions()
    centroid_x = lig_xyz[:,0].mean()
    centroid_y = lig_xyz[:,1].mean()
    centroid_z = lig_xyz[:,2].mean()
    return [centroid_x, centroid_y, centroid_z]


class LigPrep():
    # Partly adapted from https://github.com/guanjq/targetdiff/blob/main/utils/evaluation/docking_vina.py
    """Preparation of ligands.
    
    Args:
        mol (Chem.Mol): Molecule object.
        optimize (bool, optional): Generate and optimize the 3D conformation. If a conformation of the molecule is not detected, it will be forced to run. Defaults to False.
        reset_conf (bool, optional): Reset the original conformation. Defaults to False.
    """
    def __init__(self, mol:Chem.Mol, optimize=False, reset_conf=False):
        if reset_conf:
            self.mol = standard_mol(mol)
        else:
            self.mol = Chem.RemoveHs(mol) # Remove Hs to avoid error in ligprep
        self.ob_mol = pybel.readstring("mol", Chem.MolToMolBlock(self.mol))
        self.optimize = optimize
        if self.mol.GetNumConformers() == 0:
            self.optimize = True


    def ligprep(self, polaronly=False, correctforph=True, PH=7):
        """Protonate and generate a 3D conformation for the ligand.

        Args:
            polaronly (bool, optional): Add polar hydrogens only. Defaults to False.
            correctforph (bool, optional): Protonation based on pH value. Defaults to True.
            PH (int, optional): pH value. Defaults to 7.
        """
        self.ob_mol.OBMol.AddHydrogens(polaronly, correctforph, PH)
        self.mol = Chem.MolFromMolBlock(self.ob_mol.write("mol"), removeHs=False)
        if self.optimize:
            AllChem.EmbedMolecule(self.mol, AllChem.ETKDGv3())
            AllChem.MMFFOptimizeMolecule(self.mol)
        return self.mol


    def get_pdbqt(self, **kwargs):
        """Get the pdbqt string of the ligand.

        Args:
            **kwargs: Args for `LigPrep.ligprep()`.
        """
        prep_mol = self.ligprep(**kwargs)
        mk_prep = MoleculePreparation()
        molsetup_list = mk_prep(prep_mol)
        molsetup = molsetup_list[0]
        pdbqt_string = PDBQTWriterLegacy.write_string(molsetup)[0]
        return pdbqt_string


    def save_sdf(self, output_path, **kwargs):
        """Save the prepared ligand to a sdf file.

        Args:
            output_path: Path to save the sdf file.
            **kwargs: Args for `LigPrep.ligprep()`.
        """
        prep_mol = self.ligprep(**kwargs)
        write_sdf(output_path, prep_mol)
