'''
Author: Rui Qin
Date: 2025-03-15 15:56:18
LastEditTime: 2025-03-15 15:57:45
Description: 
'''
import argparse
import importlib

def main():
    parser = argparse.ArgumentParser(description="TarPass, benchmarking molecular generation tools.")
    subparsers = parser.add_subparsers(dest="module", required=True, help="Available modules")

    dock_parser = subparsers.add_parser("dock", help="Docking operations")
    dock_module = importlib.import_module("dock")
    dock_module.setup_arguments(dock_parser)

    args = parser.parse_args()
    module = importlib.import_module(args.module)
    module.execute(args)

if __name__ == "__main__":
    main()