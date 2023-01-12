# Protein Folding
A fast and simple program for folding proteins. 

## Installation
```bash
git clone https://github.com/mes3n/protein-folding
cd protein-folding
```
To use a python virtual environment:
```bash
python -m venv my_venv
source my_venv/bin/activate
```
Then install dependencies with:
```bash
pip install -r requirements.txt
```

## Dependencies
The program relies on biotite to load amino acids and chain them together (reference 1). Along with matplotlib in combinatination with biotite to graph the molecule.

## Run the Program
```bash
python bin/main.py
```
Set the number of iterations as the first argument when running the main script, default is 6:
```bash
python bin/main.py 6
```

Set the gui mode as the second argument when running the main script, default is 0 (only display molecule at the end of folding, 1 to update gui after each iteration, 2 to update gui after each fold):
```bash
python bin/main.py 6 0
```

## References
1. Kunzmann, P. & Hamacher, K. BMC Bioinformatics (2018) 19:346.
https://doi.org/10.1186/s12859-018-2367-z

2. Designing a 20-residue protein.
Neidigh, J.W., Fesinmeyer, R.M., Andersen, N.H. (2002) Nat Struct Biol 9: 425-430
