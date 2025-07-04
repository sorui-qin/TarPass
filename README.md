# Tarpass

## Setup

```Shell
conda create -n tarpass python=3.11
pip install -e . 
```

### For basic functions

```Shell
conda install -c conda-forge scipy pandas matplotlib seaborn
conda install -c conda-forge rdkit openbabel
conda install -c conda-forge biopython meeko
pip install more_itertools pyyaml tqdm networkx statsmodels scikit-posthocs

# Running PLIP
pip install plip
```

### Docking with Gnina (Suggested)

We highly suggest performancing docking with Gnina (Version 1.3).  
You can download Gnina [here](https://github.com/gnina/gnina/releases/tag/v1.3) or:

```Shell
wget https://github.com/gnina/gnina/releases/download/v1.3/gnina
```

Then add Gnina to the system path.

### Docking with AutoDock-Vina (Optional, NOT suggested)

<details>
<summary>Click to expand setup instructions</summary>

```Shell
conda install -c conda-forge swig boost-cpp libboost sphinx sphinx_rtd_theme
conda install -c conda-forge vina gemmi prody
conda install -c conda-forge autogrid # >=4.2.7 
# autogrid will conflict with umap-learn
```

### Preparation with AutoDock-Vina (Not Recommend)

```Shell
python -m pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3
```

</details>

## Input folder structure
>
> **Warning**  
> We strongly recommend that you place the molecule files to be tested according to the following folder structure. Otherwise, you will need to make the corresponding modifications in the `read_in` function of `utils/preprocess.py`.  
> We only accept SDF format files as input for 3D molecules, as reading molecules in MOL or MOL2 format may result in the loss of stereochemical information.  
> For SMILES format input, most of file formats (e.g. `.smi` or `.txt`) are acceptable).

```Text
test_folder
├── 5HT2A
│   ├── 1.sdf (or smi format, it's also accpectable to store all molecules in a single file)
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

```Text
test_folder
├── 5HT2A
│   ├── results
│   │   ├── docking_results.pkl
│   │   ├── xx_results.json
│   │   └──...
│   ├── 1.sdf
│   └── ...
└── ...
```

</details>

## Usage

```Text
tarpass [-n NUM] -p PATH {module} --option

options:
  -n NUM, --num NUM     number of unique molecules to verify per target (default: 1000).
  -p PATH, --path PATH  path to the folder where generated molecules for testing will be stored.

  {module}     modules that need to be executed: dock, interactions, dock_eval, etc.
  --option     Options specific to each module, which may vary between modules.
```

- Example

```Bash
tarpass -p test_folder dock --mode score_only
```

### Dock

**Execute docking.**

```Text
export PATH=<GNINA_PATH>:$PATH
export PYTHONUNBUFFERED=1
```

You can use the command above before you executing dock module in a slurm-based cluster.

```Text
options:
  -m , --mode {dock,score_only}  docking mode (default: dock).
  --verbose VERBOSE     verbosity level of the docking process (default: 0).
  --seed SEED           random seed for docking.
  --exhaust EXHAUST     exhaustiveness of docking.
  --poses POSES         number of poses to generate.
  --config CONFIG       path to the configuration file.
  --method METHOD       docking method to use (`gnina` or `vina`, default: 'gnina').

Optional arguments:
  --reset               reset the original 3D conformation if available
```

Default parameters for docking can be checked at `config/dock/gnina_dcok.yml`.

### DockEval

**Evaluation the docking results.**

```Text
options:
  --format {csv,json}  Output format for the evaluation results. Default in json.
```
