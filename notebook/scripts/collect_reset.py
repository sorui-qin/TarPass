import itertools, sys
from pathlib import Path
from eval.dockeval import DockEval
from utils.io import read_pkl, dump_json
from utils.preprocess import to_smiles
from utils.constant import TARGETS

def dock_eval(target_dir:Path) -> list[dict]:

    def _extract_and_eval(results) -> list[dict]:
        poses = [res['pose'] for res in results]
        scores = [res['docking score'] for res in results]
        return DockEval(poses, scores, target).evaluate()
    
    results_pkl = target_dir / 'results/gnina-reset_docking_results.pkl'
    target = target_dir.name

    if not results_pkl.exists():
        print(f"Results file {results_pkl} not found, skipping evaluation for target {target}.")
        return []
    dock_results = read_pkl(results_pkl)
    mols = [res['mol'] for res in dock_results]

    # Extract results and evaluate
    dock_eval_results = _extract_and_eval(results=dock_results)

    eval_all_results = []
    smis = to_smiles(mols)
    for idx, (smiles, dock_res) in enumerate(
        itertools.zip_longest(smis, dock_eval_results, fillvalue=None)
    ):
        info_di = {'index': idx, 'target': target, 'smiles': smiles,}
        info_di.update({'reset': dock_res})
        eval_all_results.append(info_di)
    return eval_all_results

def main(path):
    work_dir = Path(path)
    for target in TARGETS:
        print(f"Processing target: {target}")
        target_dir = work_dir / target
        results_dir = target_dir / 'results'
        eval_output = results_dir / f'reset_eval_results.json'
        if not target_dir.exists():
            print(f"Target folder {target_dir} not found, skipping evaluation.")
        dump_json(eval_output, dock_eval(target_dir))

path = sys.argv[1]
main(path)