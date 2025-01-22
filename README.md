# Computational Biology Assignment - Phylogenetic Tree Reconstruction
**Due Date: January 28, 2025**

## Overview
This code focuses on reconstructing phylogenetic trees (dendrograms) from DNA sequences, including branch lengths calculation. 

## Input Data
- Five FASTA files containing synthetic sequences of varying lengths
- Sequences generated using Jukes-Cantor rate matrix
- No indels (insertions/deletions) in the sequences
- Assumptions:
  - Position independence
  - Equal evolutionary conservation across positions
  - Time reversibility

## Task Requirements
1. Calculate evolutionary time between sequence pairs using the Jukes-Cantor model
2. Implement two tree reconstruction algorithms:
   - UPGMA (Unweighted Pair Group Method with Arithmetic Mean)
   - NJ (Neighbor Joining, Saitou-Nei method)
3. Generate phylogenetic trees in Newick format

## Output Format
- Trees should be output in Newick format
- Branch lengths should have 3 decimal places
- Include original leaf names without internal node names
- No additional output during execution
