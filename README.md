# Wildtype-and-Mutant-p53-complex-analysis
## Overview
Analyzed and modeled how the distance between the p53 and the DNA molecule changed in reference to the the wildtype complex. The significance of this project stems from sees how mutations in the structural part of the protein has effects that cascade across the whole DNA-protein complex.
<img width="3608" height="2966" alt="p53_mutation_binding_analysis" src="https://github.com/user-attachments/assets/4ed4f819-fdb8-452a-917a-927838d1276e" />


## Research Question
What are the effects of the mutations on residues 121 and 122?

## Methods
- PDB files used (3TS8.pdb, 4MZR.pdb)
- Analysis approach (distance calculations, residue contact mapping, etc.)
- Tools/languages (Python, matplotlib, numpy)

## Files
- p53mutantanalysis.py 
- 3TS8.pdb
- 4MZR.pdb
- p53_mutation_binding_analysis.png (interactive model avaliable if python file is ran)
- analysis_writeup.pdf
## Findings

DNA-contacting Residue | Wild-type | Mutant | Change
-----------------------|-----------|--------|--------
                   239 |      3.92 |   4.24 |  +0.33
                   241 |      2.69 |   2.70 |  +0.02
                   248 |      3.69 |   3.37 |  -0.32
                   273 |      2.79 |   2.91 |  +0.12
                   274 |      4.48 |   4.74 |  +0.25
                   275 |      3.48 |   3.46 |  -0.02
                   276 |      3.00 |   2.78 |  -0.22
                   277 |      3.77 |   4.11 |  +0.34
                   280 |      3.25 |   2.83 |  -0.42
                   281 |      4.62 |   4.83 |  +0.21
                   283 |      3.37 |   2.63 |  -0.74
                   284 |      4.58 |   4.43 |  -0.16
                   288 |      4.97 |   4.54 |  -0.43

## How to Run
1) Download p53mutantanalysis.py and both PDB files
2) Open the python file in an editor
3) Make sure numpy, pandas, and matplotlib are installed on your system
4) Find wildtype_file = r"C:\Users\Jason Vo\Desktop\Coding\3TS8.pdb" in line 26 and mutant_file = r"C:\Users\Jason Vo\Desktop\Coding\4MZR.pdb" in line 37
5) Replace C:\Users\Jason Vo\Desktop\Coding\3TS8.pdb with the file path of the 3TS8.pdb file in your system
6) Repeat Step 4 with the mutant file
7) Run the file
