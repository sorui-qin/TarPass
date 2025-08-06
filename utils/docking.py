'''
Author: Rui Qin
Date: 2025-03-07 19:49:34
LastEditTime: 2025-08-06 20:57:58
Description: 
'''
import numpy as np
from copy import deepcopy
from rdkit import Chem, RDLogger
from rdkit.Chem import Mol
from rdkit.Chem.rdDistGeom import EmbedMolecule
from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMolecule
from rdkit.Chem.MolStandardize.rdMolStandardize import Uncharger
RDLogger.DisableLog('rdApp.*') # type: ignore
from meeko import MoleculePreparation, PDBQTWriterLegacy
from openbabel import openbabel, pybel
openbabel.obErrorLog.SetOutputLevel(openbabel.obError)
from utils.io import write_sdf
from utils.logger import project_logger
from utils.preprocess import standard_mol

def centriod(mol:Mol) -> list[float]:
    """Calculate the centroid of a molecule.
    """
    if not mol.GetNumConformers():
        raise ValueError("The molecule has no conformation.")
    lig_xyz = mol.GetConformer().GetPositions()
    centroid_x = lig_xyz[:,0].mean()
    centroid_y = lig_xyz[:,1].mean()
    centroid_z = lig_xyz[:,2].mean()
    return [centroid_x, centroid_y, centroid_z]

def centriod_distance(mol_pred: Mol, mol_ref: Mol):
    """Calculate the distance between the centroids of two molecules."""
    centroid_pred = centriod(mol_pred)
    centroid_cond = centriod(mol_ref)
    return np.linalg.norm(
        np.array(centroid_pred) - np.array(centroid_cond)
        )

def sdf2centroid(sdf_file):
    supp = Chem.SDMolSupplier(sdf_file, sanitize=False)
    return centriod(supp[0]) if supp else None

def rdkit2obmol(mol:Mol) -> pybel.Molecule:
    return pybel.readstring("mol", Chem.MolToMolBlock(mol))

def obmol2rdkit(ob_mol:pybel.Molecule) -> Chem.Mol:
    return Chem.MolFromMolBlock(ob_mol.write("mol"), removeHs=False)

def uncharge(mol:Mol) -> Mol:
    uncharger = Uncharger()
    try:
        return uncharger.uncharge(Chem.RemoveHs(mol))
    except:
        return mol

class LigPrep():
    """Preparation of ligand.  
    The molecule will undergo hydrogen adding and protonation. 
    If the molecule has no conformation or chooses to reset the conformation, 
    conformation generation will be performed using RDKit or Open Babel and optimized using the MMFF94 force field.
    
    Args:
        mol (Chem.Mol): Molecule object.
        reset_conf (bool, optional): Reset the original conformation. Defaults to False.
    """
    def __init__(self, mol:Mol, reset_conf=False):
        self.mol = uncharge(mol) # Uncharge the molecule if formal charge != 0
        if reset_conf:
            idx = mol.GetProp('_Name')
            self.mol = standard_mol(mol)
            self.mol.SetProp('_Name', idx) # Reset the name of the molecule

    @staticmethod
    def obmol_conf(mol:Mol, minimize=False) -> pybel.Molecule:
        """Create conformation with openbabel.
        """
        ob_mol = rdkit2obmol(mol)
        ob_mol.make3D(forcefield="MMFF94", steps=50)
        if minimize:
            ob_mol.localopt(forcefield="MMFF94", steps=200)
        return ob_mol

    def ligprep(self, polaronly=True, correctforph=True, PH=7.4) -> Chem.Mol:
        """Protonate or generate a 3D conformation for the ligand.
        Args:
            polaronly (bool, optional): Add polar hydrogens only. Defaults to False.
            correctforph (bool, optional): Protonation based on pH value. Defaults to True.
            PH (int, optional): pH value. Defaults to 7.4.
        """
        p_mol = deepcopy(self.mol)

        if not self.mol.GetNumConformers(): # Check 3D conformation
            try: # RDKit pipeline
                p_mol = Chem.AddHs(p_mol, addCoords=True)
                EmbedMolecule(p_mol, useRandomCoords=True, maxAttempts=1000)
                if MMFFOptimizeMolecule(p_mol, maxIters=200): # If fail to optimize, use Open Babel
                    ob_mol = self.obmol_conf(p_mol, minimize=True)
                else:
                    ob_mol = rdkit2obmol(p_mol)
            except: # OpenBabel pipeline
                ob_mol = self.obmol_conf(p_mol, minimize=True)
        else:
            # RDmol 2 OBmol for protonation
            ob_mol = rdkit2obmol(p_mol)

        # Protonation
        backup = p_mol
        ob_mol.OBMol.AddHydrogens(polaronly, correctforph, PH)
        prep_mol = obmol2rdkit(ob_mol)
        if prep_mol is None and backup:
            project_logger.warning(f"Protonation failed in {self.mol.GetProp('_Name')}, using original conformation instead.")
            prep_mol = backup
        return prep_mol

    def get_pdbqt(self, **kwargs):
        """Get the pdbqt string of the ligand.

        Args:
            **kwargs: Args for `LigPrep.ligprep()`.
        """
        prep_mol = self.ligprep(**kwargs)
        if prep_mol:
            mk_prep = MoleculePreparation()
            molsetup_list = mk_prep(prep_mol)
            molsetup = molsetup_list[0]
            pdbqt_string = PDBQTWriterLegacy.write_string(molsetup)[0]
            return pdbqt_string
        else:
            project_logger.warning(f"Ligand preparation failed in {self.mol.GetProp('_Name')}.")
            return None

    def save_sdf(self, output_path, **kwargs):
        """Save the prepared ligand to a sdf file.

        Args:
            output_path: Path to save the sdf file.
            **kwargs: Args for `LigPrep.ligprep()`.
        """
        prep_mol = self.ligprep(**kwargs)
        if prep_mol:
            write_sdf(output_path, prep_mol)
            return True
        else:
            project_logger.warning(f"Ligand preparation failed in {self.mol.GetProp('_Name')}.")
            return False