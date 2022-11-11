# Protein Folding
A very simplified version of protein folding with help from openmm

## Installation
```bash
pip install numpy Cython
pip install molmod
```


## Installation OLD
Start by installing Anaconda using this [link](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

Make sure to add conda to path.

Setup a virtual environment with conda:
```
conda update conda
conda create -n .venv python
conda activate .venv 
```

[Install OpenMM](http://docs.openmm.org/latest/userguide/application/01_getting_started.html#installing-openmm) to the virtual environment with:
```
conda install -c conda-forge openmm 
```
If there is no need for cudatoolkit, it can be removed with:
```
conda remove --force cudatoolkit
```

## References
“T. Verstraelen, MolMod Software Library, http://molmod.ugent.be/software/”