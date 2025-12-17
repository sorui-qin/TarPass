'''
Author: Rui Qin
Date: 2025-03-15 15:56:18
LastEditTime: 2025-12-17 12:08:51
Description: 
'''
import argparse
import importlib
import sys
from pathlib import Path
from utils.io import read_yaml, temp_dir
from utils.logger import configure_third_logging

# Suppress warnings
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, message="to-Python converter for boost::shared_ptr")
warnings.filterwarnings("ignore", category=UserWarning, module="prody.utilities.misctools")
configure_third_logging()

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
    parser.add_argument('-p', '--path', required=True, type=str, help='path to the folder where generated molecules for testing will be stored.')
    parser.add_argument("-n", "--num", type=int, default=1000, help="number of unique molecules to verify per target (default: 1000).")
    subparsers = parser.add_subparsers(dest="module", required=True, help="Available modules")
    
    # Define module configurations: command -> (module_name, help_text, function_name)
    module_configs = {
        "dock": ("dock", "Docking operations", "execute"),
        "interaction": ("interaction", "Analyze interactions", "execute"),
        "eval": ("eval", "Evaluate docking and molecular properties", "eval_execute"),
        "dockeval": ("eval", "Evaluate docking results", "dockeval_execute"),
        "moleeval": ("eval", "Evaluate molecular properties", "eval_execute"),
        "analysis": ("analysis", "Analyze evalution results", "execute"),
        "collect": ("collect", "Collect evaluation results", "execute"),
    }

    parsers = {}
    for cmd, (mod_name, help_text, _) in module_configs.items():
        parsers[cmd] = subparsers.add_parser(cmd, help=help_text)

    # Lazy load module based on command line arguments
    selected_cmd = None
    for arg in sys.argv[1:]:
        if arg in module_configs:
            selected_cmd = arg
            break
            
    loaded_modules = {}
    if selected_cmd:
        mod_name = module_configs[selected_cmd][0]
        if mod_name not in loaded_modules:
            loaded_modules[mod_name] = importlib.import_module(mod_name)
        loaded_modules[mod_name].setup_arguments(parsers[selected_cmd])

    # Merge configuration for dock module
    args = parser.parse_args()
    if args.module == 'dock':
        args = merge_config(args)

    # Get the relevant module and function based on the command
    mod_name, _, function_name = module_configs[args.module]
    
    if mod_name not in loaded_modules:
        loaded_modules[mod_name] = importlib.import_module(mod_name)
    
    module = loaded_modules[mod_name]

    if hasattr(module, function_name):
        getattr(module, function_name)(args)
    else:
        raise AttributeError(f"Module {module} does not have function {function_name}")

if __name__ == "__main__":
    main()