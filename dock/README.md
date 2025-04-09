## Introduction to `Dock`

### `autobook_scripts`
The files in this directory are primarily used for the receptor preparation of the HDAC6 target and are generally not used.

### `config`
This directory contains configuration files for docking with Gnina.  
The parameters can be describle as follows:  
| Parameters | Description |
| ------ | ------ |
| Center of grid box | Center of referenced ligand |
| Size of grid box | 20Å × 20Å × 20Å |
| Number of output pose | 1 |  
Also same to AutoDock-Vina.

### `maps`
This directory contains the affinity maps for all test targets.  
We strongly recommend directly downloading and using the affinity maps from this directory when performing docking with AutoDock-Vina.  
Because affinity maps for HDAC6 should be prepared specially.  
All map files can be downloaded via [Here](todo)

### `pdbqtfiles`
This directory includes the pdbqt files for all targets.  
The pdbqt files can be used for automatical preparation, except for the affinity maps of HDAC6.