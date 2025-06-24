'''
Author: Rui Qin
Date: 2024-12-28 19:47:43
LastEditTime: 2025-06-24 21:02:39
Description: Topological and structural properties of a molecule.
'''
import copy
import networkx as nx
import numpy as np
from rdkit import Chem
from rdkit.Chem import Mol, GraphDescriptors, SpacialScore

#### Torsion and Dihedral Angles ####

def get_torsions(mol:Mol) -> list:
    # Adapted from https://github.com/fengshikun/FradNMI
    """
    Extracts the torsion angles from a molecule.
    
    Args:
        mol (Chem.rdchem.Mol): Single molecule object from RDKit.

    Returns:
        list: A list of tuples, where each tuple contains four integers representing the indices of the atoms
        that define a torsion angle in the molecule.
    """
    torsionList = []
    torsionSmarts = '[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]' # Bond between atom 2-3 should be single bond
    torsionQuery = Chem.MolFromSmarts(torsionSmarts)
    matches = mol.GetSubstructMatches(torsionQuery)
    for match in matches:
        idx2 = match[0]
        idx3 = match[1]
        bond = mol.GetBondBetweenAtoms(idx2, idx3)
        jAtom = mol.GetAtomWithIdx(idx2)
        kAtom = mol.GetAtomWithIdx(idx3)
        for b1 in jAtom.GetBonds():
            if (b1.GetIdx() == bond.GetIdx()):
                continue
            idx1 = b1.GetOtherAtomIdx(idx2)
            for b2 in kAtom.GetBonds():
                if ((b2.GetIdx() == bond.GetIdx())
                        or (b2.GetIdx() == b1.GetIdx())):
                    continue
                idx4 = b2.GetOtherAtomIdx(idx3)
                # skip 3-membered rings
                if (idx4 == idx1):
                    continue
                if mol.GetAtomWithIdx(idx4).IsInRing():
                    torsionList.append(
                        (idx4, idx3, idx2, idx1))
                    break
                else:
                    torsionList.append(
                        (idx1, idx2, idx3, idx4))
                    break
            break
    return torsionList


def get_dihedral(mol:Mol) -> list:
    # Adapted from https://github.com/gcorso/torsional-diffusion/
    """
    Extracts **all** the dihedral angles from a molecule.\n
    Compared to `get_torsion` with SMARTS matching, this function is much slower and may detect
    some *"not-so-important"* dihedral angles.

    Args:
        mol (Mol): Single molecule object from RDKit.
    Returns:
        list: A list of tuples, where each tuple contains four integers representing the indices of the atoms
        that define a dihedral angle in the molecule
    """
    dihedralList = []
    G = nx.Graph()
    for i, atom in enumerate(mol.GetAtoms()):
        G.add_node(i)
    nodes = set(G.nodes())
    for bond in mol.GetBonds():
        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        G.add_edge(start, end)
    for e in G.edges():
        G2 = copy.deepcopy(G)
        G2.remove_edge(*e)
        if nx.is_connected(G2):
            continue
        l = list(sorted(nx.connected_components(G2), key=len)[0])
        if len(l) < 2: continue
        n0 = list(G2.neighbors(e[0]))
        n1 = list(G2.neighbors(e[1]))
        dihedralList.append(
            (n0[0], e[0], e[1], n1[0])
        )
    return dihedralList


def get_torsions_number(mol:Mol, fct=get_torsions) -> int:
    return len(fct(mol))

#### Ring and Fused Ring Properties ####

class RingProp:
    """
    Calculate numerical properties of ring information.
    """
    def __init__(self, mol:Mol):
        self.info = mol.GetRingInfo()
        self.rings = self.info.AtomRings()
        # If no ring in molecule, all methods return 0
        if not self.rings:
            methods = [m for m in dir(self) if callable(getattr(self, m)) and not m.startswith("__")]
            for method in methods:
                setattr(self, method, self._return_zero)
    
    def _return_zero(self):
        return len(self.rings)

    def ring_numbers(self) -> int:
        """
        Returns:
            int: The number of rings in a molecule.
        """
        return self.info.NumRings()

    def largest_ring(self) -> int:
        """
        Returns:
            int: The number of atoms in the largest ring.
        """
        return max(len(r) for r in self.rings)
    
    def unusual(self, lower=5, upper=7, upper_only=False) -> int:
        """
        The number of *"uncommon"* rings containing either too many or too few atoms.\n
        By default, the range for the number of atoms in *"common"* rings is set between **5 and 7**,\n
        as these rings exhibit greater stability due to their ring strain.\n
        Rings that fall outside this range may be unstable and be seen as **unusual**.

        Args:
            lower (int, optional): Lower Limit of "common" rings. Defaults to 5.
            upper (int, optional): Upper Limit of "common" rings. Defaults to 7.
            upper_only (bool, optional): Only consider upper limit, accept 3 & 4-membered rings. Defaults to False.

        Returns:
            int: The number of "uncommon" rings.
        """
        if lower > 5 or upper < 7 or lower >= upper:
            raise(ValueError('Invalid input value!'))
        if lower < 3 or upper_only:
            lower = 0
        return sum([len(r) < lower or len(r) > upper for r in self.rings])


class FusedRingProp(RingProp):
    """
    Calculate numerical properties of fused ring.
    """
    def __init__(self, mol):
        super().__init__(mol)
    
    def fused_systems(self) -> list:
        """
        Returns:
            list: List with tuples of ring indexes which composing fused ring system.
        """
        n = len(self.rings)
        adj_matrix = np.zeros((n,n), dtype=int)
        for idx1 in range(n):
            for idx2 in range(idx1, n):
                if self.info.AreRingsFused(idx1, idx2):
                    adj_matrix[idx1][idx2] = 1
                    
        G = nx.from_numpy_array(adj_matrix)
        fused_sys = [g for g in nx.connected_components(G) if len(g) > 1]
        return fused_sys
    
    def largest_member(self) -> int:
        """
        Returns:
            int: The number of ring members that constitute the largest fused ring system.
        """
        return max([len(r) for r in self.fused_systems()], default=0)

#### Topological Properties ####

def bertzCT(mol:Mol) -> float:
    """
    Calculate the BertzCT (Bertz Complexity Topological) index of a molecule.
    Ref: `S. H. Bertz, J. Am. Chem. Soc., vol 103, 3599-3601 (1981).`
    """
    return GraphDescriptors.BertzCT(mol)

def spacial_score(mol:Mol) -> float:
    """
    Spacial score (SPS) is an empirical scoring system to express the spacial complexity of a compound in an uniform manner 
    and on a highly granular scale for ranking and comparison between molecules.  
    Ref: `Krzyzanowski, A.; et al. Spacial Score─A Comprehensive Topological Indicator for Small-Molecule Complexity. J. Med. Chem. 2023.`
    """
    return SpacialScore.SPS(mol)