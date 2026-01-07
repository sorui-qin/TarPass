# TarPass

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
[![ChemRxiv](https://img.shields.io/badge/ChemRxiv-Preprint-orange)](https://doi.org/10.26434/chemrxiv-2026-dhdqk)


<div align="center">
    <img src="tarpass-logo.png" style="width: 30%;">
</div>

TarPass is a comprehensive benchmark designed for target-aware *de novo* molecular generation.

## Contents

- [Quick Setup](#quick-setup)
- [Input Folder Structure](#input-folder-structure)
- [Usage](#usage)
- [CLI Help](#cli-help)
- [Notebooks and Scripts](#notebooks-and-scripts)

## Quick Setup

⚠️ We strongly recommend running TarPass on devices with a GPU.

```bash
conda env create -f tarpass.yml
conda activate tarpass
pip install -e .
```

⚠️ It should be noted that the packages required for running Jupyter Notebook are not included.

### Docking with Gnina (Suggested)

<details>
<summary>Click to expand setup instructions</summary>

We highly suggest performing docking with Gnina (Version 1.3).  
You can [download](https://github.com/gnina/gnina/releases/tag/v1.3) Gnina or:

```bash
wget https://github.com/gnina/gnina/releases/download/v1.3/gnina
```

Then add Gnina to the system path by modifying `~/.bashrc`.

```bash
export PATH=<GNINA_PATH>:$PATH
# Change <GNINA_PATH> to your own path!
export PYTHONUNBUFFERED=1
```

</details>

### Docking with AutoDock-Vina (Deprecated, NOT suggested)

<details>
<summary>Click to expand setup instructions</summary>

```bash
conda install -c conda-forge swig boost-cpp libboost sphinx sphinx_rtd_theme
conda install -c conda-forge vina gemmi prody
conda install -c conda-forge autogrid # >=4.2.7
```

### Preparation with AutoDock-Vina (Not Recommend)

```bash
python -m pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3
```

For other dependency file, please download [autodock_depends.zip](https://github.com/sorui-qin/TarPass/releases/download/v0.1.0-alpha.2/autodock_depends.zip) and copy to `dock`, then run:

```bash
cd dock
unzip autodock_depends.zip
```

</details>

## Input Folder Structure
>
> Warning  
> We strongly recommend that you place the molecule files to be tested according to the following folder structure.  
> We only accept SDF format files as input for 3D molecules, as reading molecules in MOL or MOL2 format may result in the loss of stereochemical information.
> For SMILES format input, most of file formats (e.g. `.smi` or `.txt`) are acceptable.

```text
test_folder
├── 5HT2A
│   ├── 1.sdf (or smi format,  
│   │      it's also acceptable to store all molecules in a single file)
│   ├── 2.sdf
│   └── ...
│ 
├── 5HT2A-AF
│   ├── 1.sdf
│   └── ...
└── ...
```

<details>
<summary>After the testing is completed, the file structure should be as follows:</summary>

```text
test_folder
├── 5HT2A
│   ├── results
│   │   ├── xx_docking_results.pkl
│   │   ├── xx_results.json
│   │   └──...
│   ├── 1.sdf
│   └── ...
└── ...
```

</details>

## Usage

The simplest process for benchmarking with default settings is as follows:

```text
Packaing generated molecules in designated format ->  
docking -> evaluation -> analysis -> collecting results (Optional)
```

The following sequence must be strictly followed:

### Docking

```bash
tarpass -p <path> dock
```

On devices equipped with GPUs, completing docking across all 20 proteins requires more than 16 hours, with the exact time depending on the number of molecules and the hardware performance.  

After execution, pkl files `xx-dock_docking_results.pkl` or `xx-score_only_docking_results.pkl` (3D only) containing docking scores, and docking poses will be generated under the `<target>/results` directory.

### Evaluation

```bash
tarpass -p <path> eval
```

After execution, a json file `dock_eval_results.json` with PLI-related results and a csv file `mole_eval_results.csv` for molecular properties will be generated under the `<target>/results` directory.

### Analysis

>**Warning**: [`Random_eval_results.pkl`](https://github.com/sorui-qin/TarPass/releases/download/v0.1.0-alpha/Random_eval_results.pkl) and [`Reference_eval_results.pkl`](https://github.com/sorui-qin/TarPass/releases/download/v0.1.0-alpha/Reference_eval_results.pkl) must be placed in the `data` folder.

```bash
tarpass -p <path> analysis
```

After execution, a series of csv files `analysis_<term>.csv` of analysis results will be save at `<path>`.

### Collection

```bash
tarpass -p <path> collect
```

After running, the results of eval will be collected into `xx_eval_results.pkl` under `<path>`.  

>**Note**:
>
> * The collected pickle files can also be used as input for the `analysis` module; simply replace `<path>` with the specific file path.
> * For more details on how to work with and make use of this file, please refer to `notebook/read_results.ipynb`.

### CLI Help

```text
tarpass [-n NUM] -p PATH {module} --option

options:
  -n NUM, --num NUM     number of unique molecules to verify per target (default: 1000).
  -p PATH, --path PATH  path to the folder where generated molecules for testing will be stored.

  {module}     modules that need to be executed: dock, interactions, dock_eval, etc.
  --option     Options specific to each module, which may vary between modules.
```

## Notebooks and Scripts

In addition to the main CLI tool, we provide several Jupyter Notebooks and auxiliary scripts in the `notebook/` directory for data analysis, filtering, and visualization.

### Notebooks

* `extract_nearby_res.ipynb`: A utility for extracting residues near the binding site, useful for detailed protein-ligand interaction analysis.
* `filter_with_eval.ipynb`: Provides a stepwise filtering workflow. After running the `eval` module, you can use this notebook to filter molecules based on specific targets and evaluation metrics (e.g., docking scores, PLI criteria).
* `generation_size.ipynb`: Used for analyzing the scale and distribution of the generated molecular datasets.
* `read_results.ipynb`: Instructions for using collected results.

### Auxiliary Scripts

The `notebook/scripts/` directory contains several utility scripts for specific tasks:

* `pose_check.py`: A helper script to verify and visualize docking poses.
* `cross_docking.py`: Facilitates cross-docking experiments between different targets and ligands.
* `collect_reset.py` / `tmpanaly_reset.py`: Maintenance scripts for conformation-reset docking.

> **Note**: As mentioned in the Setup section, ensure you have `jupyter` and its dependencies installed in your environment to run the `.ipynb` files, as they are not included in the default `tarpass.yml`.
