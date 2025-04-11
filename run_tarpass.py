'''
Author: Rui Qin
Date: 2025-03-15 15:56:18
LastEditTime: 2025-04-10 20:49:22
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
    temp_dir() # Clean up temp dir
    parser = argparse.ArgumentParser(description="TarPass, a target-awared molecular generation benchmarking tool.")
    subparsers = parser.add_subparsers(dest="module", required=True, help="Available modules")

    dock_parser = subparsers.add_parser("dock", help="Docking operations")
    dock_module = importlib.import_module("dock")
    dock_module.setup_arguments(dock_parser)

    args = parser.parse_args()
    module = importlib.import_module(args.module)
    args = merge_config(args)
    module.execute(args)

if __name__ == "__main__":
    main()