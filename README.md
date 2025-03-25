## Setup
```
conda create -n tarpass python=3.11
pip install -e . 
```
### For basic functions
```
conda install -c conda-forge scipy rdkit numpy openbabel biopython
pip install more_itertools

# Anlysis
conda install pandas matplotlib seaborn

# Running PILP
pip install pilp
```

### Docking with AutoDock-Vina (Optional)
```
conda install -c conda-forge swig boost-cpp libboost sphinx sphinx_rtd_theme
conda install -c conda-forge vina meeko gemmi prody
conda install -c conda-forge autogrid # >=4.2.7
```

### Preparation with AutoDock-Vina (Not Recommend)
```
python -m pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3
```

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