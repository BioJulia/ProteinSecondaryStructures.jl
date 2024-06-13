# Reference

```@meta
CollapsedDocStrings = true
```

## Basic data structure

```@docs
SSData
```

## Class retriever functions

```@docs
ss_classes
ProteinSecondaryStructures._is_class
ss_composition
ss_name
ss_code
ss_number
```

## Single-PDB runs

```@docs
stride_run
dssp_run
```

## For trajectory analysis

```@docs
ss_map
ss_content
ss_composition(::AbstractMatrix{Int}, ::Int)
```
