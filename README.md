[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://BioJulia.github.io/ProteinSecondaryStructures.jl/stable)
[![Tests](https://img.shields.io/badge/build-passing-green)](https://github.com/BioJulia/ProteinSecondaryStructures.jl/actions)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

<p align=center>
<img src=https://raw.githubusercontent.com/BioJulia/ProteinSecondaryStructures.jl/main/docs/src/assets/logo.svg>
</p>

# ProteinSecondaryStructures.jl

This package parses [STRIDE]( http://webclu.bio.wzw.tum.de/stride/) and [DSSP](https://github.com/PDB-REDO/dssp) secondary structure prediction outputs, to make them convenient to use from Julia, particularly for the analysis of MD simulations. 

## Documentation

Go to: https://BioJulia.github.io/ProteinSecondaryStructures.jl

## Citations

If you use the `STRIDE` algorithm for secondary structure prediction, please cite:

- Frishman,D & Argos,P. (1995) Knowledge-based secondary structure assignment. Proteins: structure, function and genetics, 23, 566-579.
- Kabsch,W. & Sander,C. (1983) Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features. Biopolymers, 22: 2577-2637.

If you use the `DSSP` algorithm for secondary structure prediction, please cite:

- Joosten RP, te Beek TAH, Krieger E, Hekkelman ML, Hooft RWW, Schneider R, Sander C, Vriend A series of PDB related databases for everyday needs. Nuc. Acids Res. 2010; 39:D411-D419.
- Kabsch W, Sander C. Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features. Biopolymers 1983; 22:2577-2637. 

