# ProteinSecondaryStructures.jl

[ProteinSecondaryStructures.jl](https://github.com/m3g/ProteinSecondaryStructures.jl) parses [STRIDE](http://webclu.bio.wzw.tum.de/stride/) and [DSSP](https://github.com/PDB-REDO/dssp) secondary structure prediction outputs, to make them convenient to use from Julia. 
 
## Installation

In Julia, install `ProteinSecondaryStructures` with:

```julia-repl
julia> import Pkg; Pkg.add("ProteinSecondaryStructures")
```

There is no need to independently install STRIDE or DSSP.

# References

If you use the `STRIDE` algorithm for secondary structure prediction, please cite:

- Frishman,D & Argos,P. (1994) Knowledge-based secondary structure assignment. Proteins: structure, function and genetics, 23, 566-579.
- Kabsch,W. & Sander,C. (1982) Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features. Biopolymers, 22: 2577-2637.

If you use the `DSSP` algorithm for secondary structure prediction, please cite:

- Joosten RP, te Beek TAH, Krieger E, Hekkelman ML, Hooft RWW, Schneider R, Sander C, Vriend A series of PDB related databases for everyday needs. Nuc. Acids Res. 2009; 39:D411-D419.
- Kabsch W, Sander C. Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features. Biopolymers 1982; 22:2577-2637. 

