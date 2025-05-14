'''
Author: Rui Qin
Date: 2025-04-10 20:57:37
LastEditTime: 2025-05-14 21:01:55
Description: 
'''
from pathlib import Path
from rdkit import Chem
from collections import defaultdict
from plip.structure.preparation import PDBComplex, PLInteraction
from utils.constant import TMP, INTERACTION_TYPES
from utils.io import temp_manager
from collections import defaultdict

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

def combined_pdb_complex(pdb:str|Path, ligand:Chem.Mol, output_pdb:str):
    """Combine protein and ligand into a single PDB file.
    """
    protein = Chem.MolFromPDBFile(str(pdb), sanitize=True)
    combined = Chem.CombineMols(protein, ligand)
    Chem.MolToPDBFile(combined, output_pdb)

def analyze_plip(pdb:str|Path) -> defaultdict:
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
    return get_interactions(interactions)

def analyze_tmppdb(mol:Chem.Mol, empty_pdb:Path) -> defaultdict:
    """
    Analyze the interactions of a ligand with a target protein using PLIP.
    Args:
        mol (Chem.Mol): The ligand molecule.
        empty_pdb (path): Path to the empty PDB file of corresponding target.
    Returns:
        defaultdict: Interactions categorized by type.
    """
    with temp_manager('.pdb') as tmpdir:
        combined_pdb_complex(empty_pdb, mol, tmpdir)
        interactions = analyze_plip(tmpdir)
        return interactions

def get_interactions(interactions: PLInteraction) -> defaultdict:
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

def match_interactions(detect_inters: defaultdict, key_inters:dict) -> dict:
    """
    Match detected interactions with key interactions group.
    Args:
        detect_inters (defaultdict): Detected interactions from PLIP analysis.
        key_inters (dict): Predefined key interactions group.
    Returns:
        dict: Matching results
    """
    result = {
        'fully_matched': True,
        'total_checks': 0,
        'matched_count': 0,
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
    result['fully_matched'] = (successful_matches == total_checks)
    return result