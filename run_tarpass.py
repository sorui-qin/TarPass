'''
Author: Rui Qin
Date: 2025-03-15 15:56:18
LastEditTime: 2025-06-17 17:18:57
Description: 
'''
import argparse
import importlib
from pathlib import Path
from utils.io import read_yaml, temp_dir

def merge_config(args:argparse.Namespace) -> argparse.Namespace:
    """Merge command line arguments with configuration file."""
    if args.config:
        if not Path(args.config).exists():
            raise FileNotFoundError(f"Config file '{args.config}' not found.")
        configs = read_yaml(args.config)
        if configs:
            for key, value in configs.items():
                if hasattr(args, key) and getattr(args, key, None) is None:
                    setattr(args, key, value)
    return args

def main():
    #temp_dir() # Clean up temp dir. If run in parrallel, this will cause error.
    parser = argparse.ArgumentParser(description="TarPass, a target-awared molecular generation benchmarking tool.")
    parser.add_argument("-n", "--num", type=int, default=1000, help="Number of unique molecules to check, default=1000")
    parser.add_argument('-p', '--path', required=True, type=str, help='path to the folder where generated molecules for testing will be stored.')
    subparsers = parser.add_subparsers(dest="module", required=True, help="Available modules")

    # Add subparsers for each module
    dock_parser = subparsers.add_parser("dock", help="Docking operations")
    dock_module = importlib.import_module("dock")
    dock_module.setup_arguments(dock_parser)

    interaction_parser = subparsers.add_parser("interaction", help="Analyze interactions")
    interaction_module = importlib.import_module("interaction")
    interaction_module.setup_arguments(interaction_parser)

    ### Evaluation module ###
    eval_module = importlib.import_module("eval")
    dockeval_parser = subparsers.add_parser("dockeval", help="Evaluate docking results")
    eval_module.setup_arguments(dockeval_parser)

    #moleeval_parser = subparsers.add_parser("moleeval", help="Evaluate molecular properties")
    #eval_module.setup_arguments(moleeval_parser)

    # Merge configuration for dock module
    args = parser.parse_args()
    if args.module == 'dock':
        args = merge_config(args)

    command_mapping = {
        'dock': (dock_module, 'execute'),
        'interaction': (interaction_module, 'execute'),
        'dockeval': (eval_module, 'dockeval_execute'),
        #'moleeval': (eval_module, 'mole_eval') 
    }

    # Get the relevant module and function based on the command
    module, function_name = command_mapping[args.module]

    if hasattr(module, function_name):
        getattr(module, function_name)(args)
    else:
        raise AttributeError(f"Module {module} does not have function {function_name}")

if __name__ == "__main__":
    main()