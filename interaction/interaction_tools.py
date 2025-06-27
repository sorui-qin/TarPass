'''
Author: Rui Qin
Date: 2025-04-10 20:57:37
LastEditTime: 2025-06-27 17:04:02
Description: 
'''
import json
from collections import defaultdict
from functools import partial
from pathlib import Path
from plip.structure.preparation import PDBComplex, PLInteraction
from rdkit import Chem
from rdkit.Chem import Mol
from rdkit.Chem.MolStandardize.rdMolStandardize import Uncharger
from tqdm.contrib.concurrent import process_map
from utils.constant import INTERACTION_TYPES, ROOT, TMP
from utils.io import temp_manager
from utils.preprocess import conformation_check

allkey_inters = json.load(open(ROOT/'interaction/key_interactions.json'))

def plip_tmp():
    """Create a temporary directory for PLIP analysis.
    """
    plip_tmp = TMP / 'plip'
    if not plip_tmp.exists():
        plip_tmp.mkdir(parents=True, exist_ok=True)
    else:
        for file in plip_tmp.glob('*'):
            file.unlink()
    return plip_tmp

def combined_pdb_complex(pdb:str|Path, ligand:Mol, output_pdb:str):
    """Combine protein and ligand into a single PDB file.
    """
    # Uncharge the ligand
    # Charged molecules may interfere with PLIP's judgment of hydrogen bonds
    uncharger = Uncharger()
    ligand = uncharger.uncharge(Chem.RemoveHs(ligand))
    protein = Chem.MolFromPDBFile(str(pdb), sanitize=True)
    combined = Chem.CombineMols(protein, ligand)
    Chem.MolToPDBFile(combined, output_pdb)

def analyze_plip(pdb:str|Path) -> dict[str, list]:
    """
    Execute PLIP analysis with specific pdb file.
    Args:
        pdb: Path to the PDB file.
    Returns:
        defaultdict: Interactions categorized by type.
    """
    analyzer = PDBComplex()
    analyzer.output_path = str(TMP/'plip')
    analyzer.load_pdb(str(pdb))
    analyzer.analyze()
    interactions = analyzer.interaction_sets['UNL:Z:1'] # UNL:Z:1 is default if combined with RDKit
    return dict(get_interactions(interactions))

def analyze_tmppdb(mol:Mol, empty_pdb:Path, key_inters:dict) -> dict:
    """
    Analyze the interactions of a ligand with a target protein using PLIP.
    Args:
        mol (Chem.Mol): The ligand molecule.
        empty_pdb (path): Path to the empty PDB file of corresponding target.
        key_inters (dict): Predefined key interactions group.
    Returns:
        defaultdict: Interactions categorized by type.
    """
    with temp_manager('.pdb') as tmpdir:
        combined_pdb_complex(empty_pdb, mol, tmpdir)
        detect_inters = analyze_plip(tmpdir)
    match = match_interactions(detect_inters, key_inters)
    return match

def get_interactions(interactions: PLInteraction) -> defaultdict[str, list]:
    """
    Get interactions from PLIP analysis.
    Args:
        interactions (PLInteraction): Interaction sets from PLIP analysis.
    Returns:
        defaultdict: Interactions categorized by type.
    """
    check_list = [
        interactions.hbonds_ldon,
        interactions.hbonds_pdon,
        interactions.halogen_bonds,
        interactions.hydrophobic_contacts,
        #interactions.metal_complexes, # not used
        interactions.pication_laro,
        interactions.pication_paro,
        interactions.pistacking,
        interactions.saltbridge_lneg,
        interactions.saltbridge_pneg,
        #interactions.water_bridges, # not used
    ]

    interaction_types = INTERACTION_TYPES
    inters = defaultdict(list)
    inters['H-Bond'] = []
    for name, inter in zip(interaction_types, check_list):
        for i in inter:
            inters[name].append(f'{i.reschain}_{i.restype}{i.resnr}')
    inters['H-Bond'].extend(inters['H-Bond acceptor'] + inters['H-Bond donor']) # Combine H-Bond acceptor and donor

    # PLIP excludes salt bridges from hydrogen bonds and excludes π-stacking from hydrophobic contacts,
    # so they need to be reintroduced to the corresponded group to avoid misjudgment.
    # See https://plip-tool.biotec.tu-dresden.de/plip-web/plip/help for more details
    inters['H-Bond'].extend(inters['Salt Bridge lneg'] + inters['Salt Bridge pneg'])
    inters['Hydrophobic contact'].extend(inters['Pi-stacking'])
    return inters

def match_interactions(detect_inters: dict, key_inters:dict) -> dict:
    """
    Match detected interactions with key interactions group.
    Args:
        detect_inters (defaultdict): Detected interactions from PLIP analysis.
        key_inters (dict): Predefined key interactions group.
    Returns:
        dict: Matching results
    """
    result = {
        'total_checks': 0,
        'matched_count': 0,
        'fully_matched': True,
        'matched_rate': 0,
        'detected_interactions': detect_inters,
        'matched_details': [],
        'unmatched_details': [],
    }
    total_checks, successful_matches = 0, 0

    for interaction_type, rules in key_inters.items():
        detected_resi = detect_inters.get(interaction_type, [])
        for rule in rules:
            total_checks += 1
            expected_resi = rule['residues']
            match_mode = rule['match_mode']

            if match_mode == 'all':
                is_matched = set(expected_resi).issubset(detected_resi)
            elif match_mode == 'any':
                is_matched = not set(expected_resi).isdisjoint(detected_resi)
            else:
                raise ValueError(f"Unsupported match mode: {match_mode}")

            readable_rule = {'interaction_type': interaction_type, 'rule': rule}
            if is_matched:
                successful_matches += 1
                result['matched_details'].append(readable_rule)
            else:
                result['unmatched_details'].append(readable_rule)

    result['total_checks'] = total_checks
    result['matched_count'] = successful_matches
    result['matched_rate'] = successful_matches / total_checks if total_checks > 0 else 0
    result['fully_matched'] = (successful_matches == total_checks)
    return result

def interactions(poses:list[Mol], empty_pdb:Path, key_inters:dict) -> list[dict]:
    """Anlysis the interactions of the poses with the target protein.

    Args:
        poses (list[Chem.Mol]): Mol objects with 3D conformations relative to the pocket.
        empty_pdb (Path): Path to the empty PDB file of the target protein.
        key_inters (dict): Dictionary of key interactions to analyze.

    Returns:
        list[dict]: A list of dictionaries containing the interaction analysis results for each pose.
    """
    plip_tmp() # Clean up temp dir
    if conformation_check(poses):
        # Running interaction analysis in parallel
        analyze_func = partial(analyze_tmppdb, empty_pdb=empty_pdb, key_inters=key_inters)
        return process_map(analyze_func, poses, desc="Analyzing interactions", chunksize=10)
    else:
        return []