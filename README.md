## Setup
```
conda create -n tarpass python=3.11
pip install -e . 
```
### For basic functions
```
conda install -c conda-forge scipy umap-learn
conda install -c conda-forge pandas matplotlib seaborn
conda install -c conda-forge rdkit openbabel
conda install -c conda-forge biopython meeko
pip install more_itertools pyyaml tqdm

# Running PLIP
pip install plip
```

### Docking with Gnina (suggested)
We highly suggest performancing docking with Gnina.
You can download Gnina [here](https://github.com/gnina/gnina).

### Docking with AutoDock-Vina (Optional)
<details>
<summary>Click to expand setup instructions</summary>

```
conda install -c conda-forge swig boost-cpp libboost sphinx sphinx_rtd_theme
conda install -c conda-forge vina gemmi prody
conda install -c conda-forge autogrid # >=4.2.7 
# autogrid will conflict with umap-learn
```

### Preparation with AutoDock-Vina (Not Recommend)
```
python -m pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3
```
</details>

## Input folder structure 
> **Warning**  
> We strongly recommend that you place the molecule files to be tested according to the following folder structure. Otherwise, you will need to make the corresponding modifications in the `read_in` function of `utils/preprocess.py`.  
> We only accept SDF format files as input for 3D molecules, as reading molecules in MOL or MOL2 format may result in the loss of stereochemical information.  
> For sequence format input, any file is acceptable.

```
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

After the testing is completed, the file structure should be as follows:
```
test_folder
├── 5HT2A
│   ├── results
│   │   ├── docking_results.pkl
│   │   ├── xx_results.pkl
│   │   └──...
│   ├── 1.sdf
│   └── ...
└── ...
```